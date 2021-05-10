//
// Created by ixiaohu on 2021/5/9.
//

#include <stdio.h>
#include "hfastq.h"
#include "kvec.h"
#include "utils.h"
#include "ksort.h"

typedef kvec_t(bseq1_t) bseq1_v;

typedef struct {
	ulong key;
	size_t value; // ID of the pointed read
} hash_pair_t;
typedef kvec_t(hash_pair_t) hash_pair_v;

static int usage() {
	fprintf(stderr, "Program: Reorder Reads2 by Reads1 query name\n");
	fprintf(stderr, "Usage:   reorder-R2 <reads1.fq.gz> <reads2.fq.gz> | gzip > out.fq.gz\n");
	fprintf(stderr, "         24 threads for searching by default.\n");
	fprintf(stderr, "\n");
	return 1;
}

// BKDR Hash Function
ulong BKDRHash(const char *str) {
	ulong seed = 1313; // 31 131 1313 13131 131313 etc..
	ulong hash = 0;
	while (*str) {
		hash = hash * seed + (*str++);
	}
	return hash;
}

static bseq1_v load_reads1(open_fastq_t *r1) {
	bseq1_v ret; kv_init(ret);
	while (1) {
		bseq1_t b = hfastq_fetch1(r1);
		if (b.name == NULL) break;
		int len = strlen(b.name);
		b.name[len-2] = '\0';
		free(b.comment); free(b.seq); free(b.qual);
		kv_push(bseq1_t, ret, b);
	}
	return ret;
}

static bseq1_v load_reads2(open_fastq_t *r2) {
	bseq1_v ret; kv_init(ret);
	while (1) {
		bseq1_t b = hfastq_fetch1(r2);
		if (b.name == NULL) break;
		int len = strlen(b.name);
		b.name[len-2] = '\0';
		kv_push(bseq1_t, ret, b);
	}
	return ret;
}

#define hash_pair_lt(a, b) ((a).key < (b).key)
KSORT_INIT(hash_pair, hash_pair_t, hash_pair_lt)

static hash_pair_v hash_for_reads2(const bseq1_v *r2) {
	hash_pair_v ret; kv_init(ret);
	size_t i = 0;
	for (i = 0; i < r2->n; i++) {
		const bseq1_t *p = &r2->a[i];
		hash_pair_t hp;
		hp.key = BKDRHash(p->name);
		hp.value = i;
		kv_push(hash_pair_t, ret, hp);
	}
	ks_introsort_hash_pair(ret.n, ret.a);
	size_t cnt = 0;
	for (i = 1; i < ret.n; i++) {
		if (ret.a[i].key == ret.a[i-1].key) {
			cnt++;
		}
	}
	fprintf(stderr, "The BKDRHash has %.2f %% collision rate\n", 100.0 * cnt / ret.n);
	return ret;
}

// Parallel FOR
extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

typedef struct {
	bseq1_v reads1;
	bseq1_v reads2;
	hash_pair_v hps;
} pfor_data_t;

static void search1(void *data, long seq_id, int t_id) {
	pfor_data_t *w = (pfor_data_t*) data;
	bseq1_t *b = &w->reads1.a[seq_id];
	const bseq1_v *r2 = &w->reads2;
	const hash_pair_v *hps = &w->hps;

	ulong h = BKDRHash(b->name);
	long L = 0, R = hps->n-1, mid, ans = -1;
	while (L <= R) {
		mid = (L + R) / 2;
		if (hps->a[mid].key >= h) {
			ans = mid;
			R = mid - 1;
		} else {
			L = mid + 1;
		}
	}
	long i, found = 0;
	for (i = ans; i < hps->n; i++) {
		const hash_pair_t *hp = &hps->a[i];
		if (hp->key == h) {
			const bseq1_t *p = &r2->a[hp->value];
			if (!strcmp(b->name, p->name)) {
				found = 1;
				free(b->name);
				b->name = p->name; int len = strlen(b->name);
				b->name[len] = '/'; b->name[len+1] = '2'; b->name[len+2] = '\0';
				b->comment = p->comment; b->seq = p->seq; b->qual = p->qual;
				break;
			}
		} else {
			break;
		}
	}
	if (!found) {
		fprintf(stderr, "%s didn't find the matched read\n", b->name);
		abort();
	}
}

static void output_reads2(const bseq1_v *reads1) {
	int i;
	for(i = 0; i < reads1->n; i++) {
		const bseq1_t *b = &reads1->a[i];
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
}

int reorder_R2_main(int argc, char **argv) {
	if (argc != 3) return usage();
	open_fastq_t *r1 = hfastq_open(argv[1]);
	open_fastq_t *r2 = hfastq_open(argv[2]);

	fprintf(stderr, "Stage1: Loading Reads1 query name\n");
	double realt = realtime();
	bseq1_v reads1 = load_reads1(r1);
	fprintf(stderr, "%ld sequences collected, %.2f seconds elapsed\n", reads1.n, realtime()-realt);

	fprintf(stderr, "Stage2: Loading Reads2\n");
	realt = realtime();
	bseq1_v reads2 = load_reads2(r2);
	fprintf(stderr, "%ld sequences collected, %.2f seconds elapsed\n", reads2.n, realtime()-realt);

	fprintf(stderr, "Stage3: Hash and sorting\n");
	realt = realtime();
	hash_pair_v hps = hash_for_reads2(&reads2);
	fprintf(stderr, "Hash %ld reads done, %.2f seconds elapsed\n", hps.n, realtime()-realt);

	fprintf(stderr, "Stage4: Reordering Reads2\n");
	realt = realtime();
	pfor_data_t pf_data;
	pf_data.reads1 = reads1; pf_data.reads2 = reads2; pf_data.hps = hps;
	kt_for(24, search1, &pf_data, reads1.n);
	fprintf(stderr, "Reordering done, %.2f seconds elapsed\n", realtime()-realt);

	fprintf(stderr, "Stage5: Outputting reordered Reads2\n");
	realt = realtime();
	output_reads2(&reads1); // Note that the reordered reads2 saved in reads1.
	fprintf(stderr, "Done, %.2f seconds elapsed\n", realtime()-realt);

	free(reads1.a); free(reads2.a); free(hps.a);

	hfastq_close(r1);
	hfastq_close(r2);
	return 0;
}