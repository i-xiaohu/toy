//
// Created by ixiaohu on 2021/3/7.
//

#include <stdio.h>
#include <stdint.h>

#include "hfastq.h"
#include "kvec.h"

static int usage() {
	fprintf(stderr, "Program: sort_reads <in.fq> <reads> <out.fq>\n");
	return 1;
}

typedef kvec_t(int) int_v;

typedef struct {
	int id; // -1 means not end of a read.
	int_v more_ids; // todo: 如果有更多匹配项，遍历之后请释放这里的内存
	long sons[5];
} trie_node_t;
typedef kvec_t(trie_node_t) trie_node_v;
trie_node_v trie_nodes;

typedef kvec_t(bseq1_t) read_v;
read_v reads;

static trie_node_t new_node(int id) {
	trie_node_t ret;
	ret.id = id;
	kv_init(ret.more_ids);
	memset(ret.sons, -1, sizeof(ret.sons));
	return ret;
}

static void insert_string(trie_node_v *nodes, int id, int len, const char *s) {
	int i;
	long p = 0;
	fprintf(stderr, "%s %d\n", s, len);
	for(i = 0; i < len; i++) {
		int c = 4;
		if(s[i] == 'A') c = 0;
		else if(s[i] == 'C') c = 1;
		else if(s[i] == 'G') c = 2;
		else if(s[i] == 'T') c = 3;
		long child = nodes->a[p].sons[c];
		if(child == -1) {
			child = nodes->n;
			nodes->a[p].sons[c] = child;
			if(i != len-1) kv_push(trie_node_t, *nodes, new_node(-1));
			else kv_push(trie_node_t, *nodes, new_node(id));
		} else {
			if(i == len-1) {
				if(nodes->a[child].id == -1) nodes->a[child].id = id;
				else kv_push(int, nodes->a[child].more_ids, id);
			}
		}
		fprintf(stderr, "[%ld] -> [%ld] [%d] [%d, %c]\n", p, child, nodes->a[child].id, c, s[i]);
		p = child;
	}
}

static void dump_sorted_reads(long fa) {
	const trie_node_t *p = &trie_nodes.a[fa];
	int i;
	if(p->id != -1) {
		printf("%s\n", reads.a[p->id].seq);
		for(i = 0; i < p->more_ids.n; i++) {
			printf("%s\n", reads.a[p->more_ids.a[i]].seq);
		}
	}
	for(i = 0; i < 5; i++) {
		long child = p->sons[i];
		if(child != -1) {
			dump_sorted_reads(child);
		}
	}
}

int sort_reads_main(int argc, char *argv[]) {
	if(argc != 4) return usage();
	open_fastq_t *fi = hfastq_open(argv[1]);
	if(fi == NULL) {
		fprintf(stderr, "Open %s files\n", argv[1]);
		return 1;
	}

	kv_init(reads);
	kv_init(trie_nodes);
	kv_push(trie_node_t, trie_nodes, new_node(-1)); // Root node
	while(1) {
		bseq1_t read = hfastq_fetch1(fi);
		if(read.l_seq == 0) break;
		insert_string(&trie_nodes, reads.n, read.l_seq, read.seq);
		kv_push(bseq1_t, reads, read);
	}
	dump_sorted_reads(0);

	return 0;
}