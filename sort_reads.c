//
// Created by ixiaohu on 2021/3/7.
//

#include <stdio.h>
#include <stdint.h>

#include "hfastq.h"
#include "kvec.h"
#include "ksort.h"

static int usage() {
	fprintf(stderr, "Program: sort_reads <in.fq> <reads> <out.fq>\n");
	return 1;
}

typedef struct {
	bseq1_t seq;
	uint64_t hash_v[5];
} read_t;

typedef kvec_t(read_t) read_v;
#define read_lt(a, b) (((a).hash_v[0]) < ((b).hash_v[0]) \
	|| ((a).hash_v[1]) < ((b).hash_v[1])                          \
	|| ((a).hash_v[2]) < ((b).hash_v[2])                          \
	|| ((a).hash_v[3]) < ((b).hash_v[3])                          \
	|| ((a).hash_v[4]) < ((b).hash_v[4])                          \
	|| (strcmp((a).seq.seq, (b).seq.seq) < 0))
KSORT_INIT(read, read_t, read_lt)

static uint64_t base2num(char x) {
	if(x == 'A') return 0;
	else if(x == 'C') return 1;
	else if(x == 'G') return 2;
	else if(x == 'T') return 3;
	else return 0;
}

static void set_hashv(int len, const char *a, uint64_t v[]) {
	int i, j = 0;
	// The len should less 160 stored five 64-bit number
	for(i = 0; i < len; i++) {
		v[j] += base2num(a[i]);
		v[j] <<= 2U;
		if((i+1) % 32 == 0) j++;
	}
}

int sort_reads_main(int argc, char *argv[]) {
	if(argc != 4) return usage();
	open_fastq_t *fi = hfastq_open(argv[1]);
	if(fi == NULL) {
		fprintf(stderr, "Open %s files\n", argv[1]);
		return 1;
	}

	read_v reads; kv_init(reads);
	read_t tmp;
	while(1) {
		tmp.seq = hfastq_fetch1(fi);
		if(tmp.seq.l_seq == 0) break;
		memset(tmp.hash_v, 0, sizeof(tmp.hash_v));
		set_hashv(tmp.seq.l_seq, tmp.seq.seq, tmp.hash_v);
		kv_push(read_t, reads, tmp);
	}

	ks_introsort(read, reads.n, reads.a);
	int i, cnt = 0;
	for(i = 1; i < reads.n; i++) {
		if(strcmp(reads.a[i-1].seq.seq, reads.a[i].seq.seq) <= 0) {
			cnt++;
		}
	}
	printf("[%d, %ld]", cnt+1, reads.n);
	return 0;
}