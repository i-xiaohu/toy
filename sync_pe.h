//
// Created by 63175 on 2019/11/22.
//

#ifndef TOY_1_1_SYNC_PE_H
#define TOY_1_1_SYNC_PE_H

typedef struct {
	int l_seq, id;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

#ifdef __cplusplus
extern "C" {
#endif

	bseq1_t *bseq_read1(int chunk_size, int *n_, long *bytes, void *ks1);

#ifdef __cplusplus
}
#endif

#endif //TOY_1_1_SYNC_PE_H
