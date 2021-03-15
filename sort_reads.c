//
// Created by ixiaohu on 2021/3/7.
//

#include <stdio.h>

#include "hfastq.h"
#include "kvec.h"
#include "utils.h"
#include "kstring.h"

#define TRIE_CHAR 1000000000 // At most 1G characters constructing a trie.

static int usage() {
	fprintf(stderr, "Program: sort_reads <source.fq> <reordered reads>\n");
	fprintf(stderr, "Usage: sort-reads in.fq.gz in.reads.gz | gzip > out.fq.gz\n");
	return 1;
}

typedef kvec_t(int) int_v;

typedef struct {
	int id; // -1 means not end of a read.
	int_v more_ids;
	long sons[5];
} trie_node_t;
typedef kvec_t(trie_node_t) trie_node_v;
trie_node_v trie_nodes;

typedef kvec_t(bseq1_t) read_v;
read_v source_reads, reorder_reads;

/** Return formatted numbers and time.
 * 12345345 -> 12,345,345
 * 12898 -> 03:34:58 */
typedef kvec_t(kstring_t) kstring_v;
kstring_v fmt_nums;

static char* fmtn(long num) {
	kstring_t *ret = kv_pushp(kstring_t, fmt_nums);
	memset(ret, 0, sizeof(*ret));
	int i, cnt = 0, neg = 0;
	if(num < 0) neg = 1, num = -num;
	do {
		kputc('0'+(num%10), ret);
		cnt++;
		num /= 10;
		if(cnt == 3 && num != 0) {
			kputc(',', ret);
			cnt = 0;
		}
	} while(num);
	if(neg) kputc('-', ret);
	for(i = 0; i < ret->l/2; i++) {
		char c = ret->s[ret->l-1-i];
		ret->s[ret->l-1-i] = ret->s[i];
		ret->s[i] = c;
	}
	return ret->s;
}

static char* fmtt(int sec) {
	kstring_t *ret = kv_pushp(kstring_t, fmt_nums);
	memset(ret, 0, sizeof(*ret));
	int hour = sec / 3600; sec -= hour * 3600;
	int min = sec / 60; sec -= min * 60;
	ksprintf(ret,"%02d:%02d:%02d", hour, min, sec);
	return ret->s;
}

static void fn_destroy() {
	int i;
	for(i = 0; i < fmt_nums.n; i++) {
		kstring_t *p = &fmt_nums.a[i];
		free(p->s);
	}
	free(fmt_nums.a);
}

/** Construct trie. */
static trie_node_t new_node(int id) {
	trie_node_t ret;
	ret.id = id;
	kv_init(ret.more_ids);
	memset(ret.sons, -1, sizeof(ret.sons));
	return ret;
}

static inline int char2num(char x) {
	int c;
	if(x == 'A') c = 0;
	else if(x == 'C') c = 1;
	else if(x == 'G') c = 2;
	else if(x == 'N') c = 3;
	else if(x == 'T') c = 4;
	else {
		fprintf(stderr, "[%s] Do not recognize the character '%c'\n", __func__, x);
		abort();
	}
	return c;
}

int same_cnt;
static void insert_string(trie_node_v *nodes, int id, int len, const char *s) {
	int i;
	long p = 0;
	for(i = 0; i < len; i++) {
		int c = char2num(s[i]);
		long child = nodes->a[p].sons[c];
		if(child == -1) {
			child = nodes->n;
			nodes->a[p].sons[c] = child;
			if(i != len-1) kv_push(trie_node_t, *nodes, new_node(-1));
			else kv_push(trie_node_t, *nodes, new_node(id));
		} else {
			if(i == len-1) {
				if(nodes->a[child].id == -1) nodes->a[child].id = id;
				else {
					same_cnt++;
					kv_push(int, nodes->a[child].more_ids, id);
				}
			}
		}
		p = child;
	}
}

static void destroy_trie(trie_node_v *nodes) {
	int i;
	for(i = 1; i < nodes->n; i++) {
		const trie_node_t *p = &nodes->a[i];
		if(p->more_ids.a) free(p->more_ids.a);
	}
	free(nodes->a);
}

static void fill_read(bseq1_t *dst, const bseq1_t* src) {
	dst->name = src->name;
	dst->comment = src->comment;
	dst->qual = src->qual;
}

static int search(const trie_node_v *nodes, bseq1_t *r) {
	int i;
	long p = 0;
	for(i = 0; i < r->l_seq; i++) {
		int c = char2num(r->seq[i]);
		const trie_node_t *fa = &nodes->a[p];
		if(fa->sons[c] == -1) return -1; // Not matched.
		trie_node_t *son = &nodes->a[fa->sons[c]];
		if(i == r->l_seq - 1) {
			// Consider ids as stack.
			if(son->id == -1) return -1; // Not end-to-end matched.
			if(son->id == -2) return -2; // Matched at end, but token by previous reads.
			if(son->more_ids.n > 0) {
				int id = son->more_ids.a[--son->more_ids.n];
				fill_read(r, &source_reads.a[id]);
			} else {
				int id = son->id;
				fill_read(r, &source_reads.a[id]);
				son->id = -2;
			}
		}
		p = fa->sons[c];
	}
	return 0;
}

int sort_reads_main(int argc, char *argv[]) {
	if(argc != 3) return usage();
	double rtime = realtime();
	kv_init(fmt_nums);
	open_fastq_t *sf = hfastq_open(argv[1]);
	if(sf == NULL) {
		fprintf(stderr, "Open source FASTQ file %s failed\n", argv[1]);
		return 1;
	}
	open_fastq_t *rr = hfastq_open(argv[2]);
	if(rr == NULL) {
		fprintf(stderr, "Open reordered reads file %s failed\n", argv[2]);
		return 1;
	}

	kv_init(reorder_reads);
	long bases_memory = 0;
	while(1) {
		bseq1_t read = hfastq_fetch1(rr);
		if(read.l_seq == 0) break;
		free(read.name); free(read.comment);
		read.name = NULL;
		bases_memory += read.l_seq;
		kv_push(bseq1_t , reorder_reads, read);
	}
	fprintf(stderr, "Input %s reordered reads[%sMB], %s elapsed.\n",
		 fmtn(reorder_reads.n), fmtn(bases_memory/1024/1024), fmtt((int)realtime()-rtime));

	// Constructing trie within chunk size.
	fprintf(stderr, "Starting reading source FASTQ file '%s'\n", argv[1]);
	int i, done = 0, chunk_n = 1, all_matched = 0;
	same_cnt = 0;
	kv_init(source_reads); // For keeping qname and qual.
	while(!done) {
		kv_init(trie_nodes);
		kv_push(trie_node_t, trie_nodes, new_node(-1)); // Root node
		long reads_memory = 0, bases_n = 0, reads_n = 0;
		while(1) {
			bseq1_t read = hfastq_fetch1(sf);
			if(read.l_seq == 0) {
				done = 1;
				break;
			}
			reads_n++;
			insert_string(&trie_nodes, source_reads.n, read.l_seq, read.seq);
			kv_push(bseq1_t, source_reads, read);
			free(read.seq); // The bases sequence is useless after inserted into trie.
			reads_memory += strlen(read.name);
			if(read.comment) reads_memory += strlen(read.comment);
			bases_n += read.l_seq;
			reads_memory += read.l_seq; // Quality score.
			reads_memory += sizeof(read);
			if(bases_n >= TRIE_CHAR) break;
		}
		fprintf(stderr, "Chunk %d: constructed trie[%sMB] of by %s reads[%sMB], %s elapsed.\n",
		  chunk_n++,
		  fmtn(bases_n * sizeof(trie_node_t)/1024/1024),
		  fmtn(reads_n),
		  fmtn(reads_memory/1024/1024),
		  fmtt((int)realtime()-rtime));

		// Search
		int matched = 0;
		for(i = 0; i < reorder_reads.n; i++) {
			bseq1_t *p = &reorder_reads.a[i];
			if(p->name) continue; // Found the matched in trie.
			if(search(&trie_nodes, p) == 0) {
				matched++;
			}
		}
		all_matched += matched;
		fprintf(stderr, "%s read matched.\n", fmtn(matched));

		// Destroy trie.
		destroy_trie(&trie_nodes);
		if(matched != reads_n) fprintf(stderr, "Warning: The trie is not matched completely.\n");
	}

	fprintf(stderr, "Found %s same reads, taking %.2f%% of the reads in trie.\n", fmtn(same_cnt), 100.0*same_cnt/source_reads.n);
	if(all_matched != reorder_reads.n) {
		fprintf(stderr, "Warning: %d reordered reads are not matched.\n", reorder_reads.n - all_matched);
	} else {
		fprintf(stderr, "All of reordered reads are matched.\n");
	}

	// Output complete reordered FASTQ.
	for(i = 0; i < reorder_reads.n; i++) {
		const bseq1_t *b = &reorder_reads.a[i];
		if(b->qual) fprintf(stdout, "@%s", b->name);
		else fprintf(stdout, ">%s", b->name);
		if(b->comment) fprintf(stdout, " %s", b->comment);
		fprintf(stdout, "\n");
		fprintf(stdout, "%s\n", b->seq);
		if(b->qual) {
			fprintf(stdout, "+\n");
			fprintf(stdout, "%s\n", b->qual);
		} // else FASTA.
	}

	// Free allocated memory.
	free(source_reads.a);
	for(i = 0; i < reorder_reads.n; i++) {
		const bseq1_t *p = &reorder_reads.a[i];
		free(p->name);free(p->comment); free(p->seq); free(p->qual);
	}
	free(reorder_reads.a);
	hfastq_close(sf); hfastq_close(rr);
	fn_destroy();

	return 0;
}