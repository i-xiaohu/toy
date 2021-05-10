//
// Created by ixiaohu on 2021/5/9.
//

#include <stdio.h>
#include "hfastq.h"
#include "kvec.h"
#include "utils.h"
#include "ksort.h"

typedef kvec_t(bseq1_t) bseq1_v;

static int usage() {
	fprintf(stderr, "Program: Reorder Reads2 by Reads1 query name\n");
	fprintf(stderr, "Usage:   reorder-R2 <reads1.fq.gz> <reads2.fq.gz> | gzip > out.fq.gz\n");
	fprintf(stderr, "         24 threads for searching by default.\n");
	fprintf(stderr, "\n");
	return 1;
}


// BKDR Hash Function
int BKDRHash(const char *str) {
	unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
	unsigned int hash = 0;
	while (*str) {
		hash = hash * seed + (*str++);
	}
	return (hash & 0x7FFFFFFF);
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

#define bseq_lt(a, b) ((a).id < (b).id)
KSORT_INIT(bseq, bseq1_t, bseq_lt)

static bseq1_v load_reads2(open_fastq_t *r2) {
	bseq1_v ret; kv_init(ret);
	while (1) {
		bseq1_t b = hfastq_fetch1(r2);
		if (b.name == NULL) break;
		int len = strlen(b.name);
		b.name[len-2] = '\0';
		b.id = BKDRHash(b.name);
		kv_push(bseq1_t, ret, b);
	}
	ks_introsort_bseq(ret.n, ret.a);
	int i, cnt = 0;
	for (i = 1; i < ret.n; i++) {
		if (ret.a[i].id == ret.a[i-1].id) {
			cnt++;
		}
	}
	fprintf(stderr, "The BKDR hash has %.2f collisions\n", 100.0 * cnt / ret.n);
	return ret;
}

// Parallel FOR
extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

typedef struct {
	bseq1_v reads1;
	bseq1_v reads2;
} pfor_data_t;

static void search1(void *data, long seq_id, int t_id) {
	pfor_data_t *w = (pfor_data_t*) data;
	bseq1_t *b = &w->reads1.a[seq_id];
	const bseq1_v *r2 = &w->reads2;

	int h = BKDRHash(b->name);
	int L = 0, R = r2->n-1, mid, ans = -1;
	while (L <= R) {
		mid = (L + R) / 2;
		if (r2->a[mid].id >= h) {
			ans = mid;
			R = mid - 1;
		} else {
			L = mid + 1;
		}
	}
	int i, found = 0;
	for (i = ans; i < r2->n; i++) {
		const bseq1_t *p = &r2->a[i];
		if (p->id == h) {
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
	fprintf(stderr, "%.2f seconds elapsed\n", realtime()-realt);

	fprintf(stderr, "Stage2: Loading, hash and sort Reads2\n");
	realt = realtime();
	bseq1_v reads2 = load_reads2(r2);
	fprintf(stderr, "%.2f seconds elapsed\n", realtime()-realt);

	fprintf(stderr, "Stage3: Reordering Reads2\n");
	realt = realtime();
	pfor_data_t pf_data;
	pf_data.reads1 = reads1;
	pf_data.reads2 = reads2;
	kt_for(24, search1, &pf_data, reads1.n);
	fprintf(stderr, "%.2f seconds elapsed\n", realtime()-realt);

	fprintf(stderr, "Stage4: Outputting reordered Reads2\n");
	realt = realtime();
	output_reads2(&reads1); // Note that the reordered reads2 saved in reads1.
	fprintf(stderr, "%.2f seconds elapsed\n", realtime()-realt);

	free(reads1.a); free(reads2.a);

	hfastq_close(r1);
	hfastq_close(r2);
	return 0;
}