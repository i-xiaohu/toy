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
	fprintf(stderr, "Program: sort_reads [options] <source.fq> <reordered reads>\n");
	fprintf(stderr, "Usage:   sort-reads in.fq.gz reordered.reads | gzip > out.fq.gz\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -t [INT] Thread number to search[24]\n");
	fprintf(stderr, "\n");
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

/** Parallel supported */
extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

typedef kvec_t(long) long_v;

typedef struct {
	long_v *sub_rid, *sub_tid;  // read_id and tree_node_id.
	long_v rid, tid;
} worker_t;

static void search(void *data, long seq_id, int t_id) {
	worker_t *w = (worker_t*)data;
	bseq1_t *r = &reorder_reads.a[seq_id];
	long_v *rid = &w->sub_rid[t_id], *tid = &w->sub_tid[t_id];
	if(r->name) return ;
	int i;
	long p = 0;
	for(i = 0; i < r->l_seq; i++) {
		int c = char2num(r->seq[i]);
		const trie_node_t *fa = &trie_nodes.a[p];
		if(fa->sons[c] == -1) return ; // The query base has not matched node.
		trie_node_t *son = &trie_nodes.a[fa->sons[c]];
		if(i == r->l_seq - 1) {
			if(son->id == -1) return ; // The ending base is matching 'son', but 'son' is not the end of any read.
			kv_push(long, *rid, seq_id);
			kv_push(long, *tid, fa->sons[c]);
		}
		p = fa->sons[c];
	}
}

static int sumup_search(int n_threads, worker_t *w) {
	int i, j;
	kv_init(w->rid); kv_init(w->tid);
	for(i = 0; i < n_threads; i++) {
		const long_v *sr = &w->sub_rid[i];
		const long_v *st = &w->sub_tid[i];
		for(j = 0; j < sr->n; j++) {
			kv_push(long, w->rid, sr->a[j]);
			kv_push(long, w->tid, st->a[j]);
		}
		free(sr->a); free(st->a);
	}
	free(w->sub_rid); free(w->sub_tid);

	int ret = 0;
	for(i = 0; i < w->rid.n; i++) {
		long seq_id = w->rid.a[i];
		long son_id = w->tid.a[i];
		bseq1_t *r = &reorder_reads.a[seq_id];
		trie_node_t *son = &trie_nodes.a[son_id];
		if(son->id == -2) {
			continue; // The stack is empty, the read has no match.
		}
		// See the ids as a stack.
		if(son->more_ids.n > 0) {
			int id = son->more_ids.a[--son->more_ids.n];
			fill_read(r, &source_reads.a[id]);
		} else {
			int id = son->id;
			fill_read(r, &source_reads.a[id]);
			son->id = -2; // Set the stack empty.
		}
		ret++;
	}
	free(w->rid.a); free(w->tid.a);
	return ret;
}


int sort_reads_main(int argc, char *argv[]) {
	int c, n_threads = 24;
	while((c = getopt(argc, argv, "t:")) >= 0) {
		if(c == 't') n_threads = atoi(optarg);
		else return usage();
	}
	if(argc-optind != 2) return usage();

	kv_init(fmt_nums);

	double rtime = realtime();
	fprintf(stderr, "Loading reordered reads files '%s' into memory.\n", argv[2]);
	gzFile rr = gzopen(argv[2], "r");
	if(rr == NULL) {
		fprintf(stderr, "Open reordered reads file %s failed\n", argv[2]);
		return 1;
	}
	kv_init(reorder_reads);
	long bases_memory = 0;
	char buf[10240];
	while(gzgets(rr, buf, sizeof(buf))) {
		int len = strlen(buf);
		if(buf[len-1] == '\n') buf[--len] = '\0';
		bseq1_t read; memset(&read, 0, sizeof(read));
		read.seq = strdup(buf);
		read.l_seq = len;
		bases_memory += read.l_seq;
		kv_push(bseq1_t , reorder_reads, read);
	}
	fprintf(stderr, "Input %s reordered reads[%sMB], %s elapsed.\n",
		 fmtn(reorder_reads.n), fmtn(bases_memory/1024/1024), fmtt((int)realtime()-rtime));

	// Constructing trie within chunk size.
	fprintf(stderr, "Starting reading source FASTQ file '%s'\n", argv[1]);
	open_fastq_t *sf = hfastq_open(argv[1]);
	if(sf == NULL) {
		fprintf(stderr, "Open source FASTQ file %s failed\n", argv[1]);
		return 1;
	}
	int i, done = 0, chunk_n = 1, all_matched = 0;
	same_cnt = 0;
	kv_init(source_reads); // For keeping qname and qual.
	while(!done) {
		// Construct trie.
		double time1 = realtime();
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
		  fmtt((int)realtime()-time1));

		// Search
		time1 = realtime();
		worker_t w;
		w.sub_tid = calloc(n_threads, sizeof(long_v));
		w.sub_rid = calloc(n_threads, sizeof(long_v));
		kt_for(n_threads, search, &w, reorder_reads.n);
		int matched = sumup_search(n_threads, &w);
		all_matched += matched;
		fprintf(stderr, "%s read matched, %s elapsed.\n", fmtn(matched), fmtt((int)realtime()-time1));
		if(matched != reads_n) fprintf(stderr, "Warning: The trie is not matched completely.\n");

		// Destroy trie.
		destroy_trie(&trie_nodes);
	}

	fprintf(stderr, "Found %s same reads, taking %.2f%% of the reads in trie.\n", fmtn(same_cnt), 100.0*same_cnt/source_reads.n);
	if(all_matched != reorder_reads.n) {
		fprintf(stderr, "Warning: %d reordered reads are not matched.\n", reorder_reads.n - all_matched);
	} else {
		fprintf(stderr, "All of reordered reads are matched.\n");
	}

	// Output complete reordered FASTQ and free allocated memory.
	free(source_reads.a);
	for(i = 0; i < reorder_reads.n; i++) {
		const bseq1_t *b = &reorder_reads.a[i];
		if(b->qual) fprintf(stdout, "@%s", b->name); // FASTQ
		else fprintf(stdout, ">%s", b->name); // FASTA.
		if(b->comment) fprintf(stdout, " %s", b->comment);
		fprintf(stdout, "\n");
		fprintf(stdout, "%s\n", b->seq);
		if(b->qual) {
			fprintf(stdout, "+\n");
			fprintf(stdout, "%s\n", b->qual);
		}
		free(b->name); free(b->comment); free(b->seq); free(b->qual);
	}
	free(reorder_reads.a);
	hfastq_close(sf); gzclose(rr);
	fn_destroy();

	return 0;
}