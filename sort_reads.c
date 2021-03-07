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
#define read_lt(a, b) (((a).hash_v) < ((b).hash_v) || (strcmp((a).seq.seq, (b).seq.seq) < 0))
KSORT_INIT(read, read_t, read_lt)

static uint64_t base2num(char x) {
	if(x == 'A') return 0;
	else if(x == 'C') return 1;
	else if(x == 'G') return 2;
	else if(x == 'T') return 3;
	else return 0;
}

static void set_hashv(int len, const char *a, uint64_t v[]) {
	int i; uint64_t ret = 1;
	// The len should less 160 stored five 64-bit number
	for(i = 0; i < len; i++) {
		ret += base2num(a[i]);
		ret <<= 2U;
	}
}

typedef kvec_t(uint64_t) uint64_v;
#define uint64_lt(a, b) ((a) < (b))
KSORT_INIT(uint64, uint64_t, uint64_lt)

int sort_reads_main(int argc, char *argv[]) {
	if(argc != 4) return usage();
	open_fastq_t *fi = hfastq_open(argv[1]);
	if(fi == NULL) {
		fprintf(stderr, "Open %s files\n", argv[1]);
		return 1;
	}
	read_t tmp;

	uint64_v p; kv_init(p);
	do {
		tmp.seq = hfastq_fetch1(fi);
		memset(tmp.hash_v, 0, sizeof(tmp.hash_v));
		set_hashv(tmp.seq.l_seq, tmp.seq.seq, tmp.hash_v);
		kv_push(uint64_t, p, tmp.hash_v[0]);
	} while(tmp.seq.l_seq > 0);
	ks_introsort(uint64, p.n, p.a);
	int i, hit = 0;
	for(i = 1; i < p.n; i++) {
		if(p.a[i] == p.a[i-1]) {
			hit++;
		}
	}
	fprintf(stderr, "Hit_ratio[%%] %.2f\n", 100.0 * hit / p.n);

	return 0;
}