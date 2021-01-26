//
// Created by ixiaohu on 2021/1/24.
//

#ifndef TOY_HFASTQ_H
#define TOY_HFASTQ_H

#include "zlib.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int fd;
	gzFile fp;
	void *ko;
	kseq_t *ks;
} open_fastq_t;

typedef struct {
	int l_seq, id;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

open_fastq_t* hfastq_open(const char *fn);
void hfastq_close(open_fastq_t *p);
bseq1_t hfastq_fetch1(kseq_t *ks);

#endif //TOY_HFASTQ_H
