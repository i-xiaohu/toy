//
// Created by 63175 on 2019/8/1.
//

#include <stdio.h>
#include <unistd.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include "samop.h"

void sam2sfqUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    sam2sfq get sorted FASTQ file from SAM file\n");
	fprintf(stderr, "Usage:      sam2sfq [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -i STR        SAM format input file\n");
	fprintf(stderr, "            -o STR        SFQ format output file\n");
	fprintf(stderr, "\n");
}

int doSAMToSFQ(FILE *fi, FILE *fo) {
	SAMInfo_t info = samGetInfo(fi);
	SAMNode_v *sams = &info.rec1;
	SAMOpAux_t *aux = &info.aux;
	kstring_t *bs = &aux->bigStr;
	char *SO = samGetHeaderValue(&info, "SO");
	if(SO==NULL || strcmp(SO, "coordinate") != 0) {
		fprintf(stderr, "SAM文件没有按坐标排序\n");
		return 1;
	}
	// FIXME: 去冗余的策略不好，chimeric reads不能处理

	/* 开始输出sfq文件，已经完成去冗余工作 */
	int  INF=1000000000; // 染色体长度不超过20,000,000
	size_t i, j;
	for(i=0; i<sams->n; ++i) {
		SAMNode_t *sam = &sams->a[i];
		/* 不需要维持原fastq中元数据的状态 */
		fprintf(fo, "@%s\n", bs->s+sam->qname); // 需要加上@做开头

		/* SAM中的SEQ. FIXME：这里的序列可能被hard clipped, 序列可能不全 */
		char *seq = bs->s + sam->seq;
		size_t sLen = strlen(seq);
		for(j=0; j<sLen; ++j) {
			if(seq[j] == 'N') {
				seq[j] = 'A';
			}
		}
		fprintf(fo, "%s\n", bs->s+sam->seq);

		/* 输出两条序列之间的距离 */
		if(i == 0) {
			fprintf(fo, "%d\n", INF);
		} else {
			SAMNode_t *pre = &sams->a[i-1];
			if(strcmp(bs->s+pre->rname, bs->s+sam->rname) == 0) {
				int dis = sam->pos - pre->pos;
				fprintf(fo, "%d\n", dis); // 先不管两条序列之间的错误数，只要相邻即可。
			} else {
				fprintf(fo, "%d\n", INF);
			}
		}

		/* 使出质量分数 */
		fprintf(fo, "%s\n", bs->s+sam->qual);
	}
	if(fi != NULL) fclose(fi);
	if(fo != NULL) fclose(fo);
	return 0;
}

int sam2sfq_main(int argc, char *argv[]) {
	if(argc == 1) {
		sam2sfqUsage();
		return 0;
	}
	int c;
	FILE *fSAM = NULL, *fSFQ = NULL;
	while((c=getopt(argc, argv, "i:o:")) > -1) {
		if(c == 'i') {
			fSAM = fopen(optarg, "r");
			if(fSAM == NULL) {
				fprintf(stderr, "fail to open SAM file\n");
				return 1;
			}
		} else if(c == 'o') {
			fSFQ = fopen(optarg, "w");
			if(fSFQ == NULL) {
				fprintf(stderr, "fail to open SFQ file\n");
				return 1;
			}
		} else {
			fprintf(stderr, "option error!\n");
			sam2sfqUsage();
			return 1;
		}
	}
	return doSAMToSFQ(fSAM, fSFQ);
}