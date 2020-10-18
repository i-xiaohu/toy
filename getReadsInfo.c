//
// Created by 63175 on 2019/7/31.
//

/*
 * 获得序列相关的信息
 */

#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

void getReadsInfoUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    getReadsInfo get reads information of FASTQ file\n");
	fprintf(stderr, "Usage:      SAMCheck [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -f STR        FASTQ file name\n");
	fprintf(stderr, "\n");
}

void getReadsInfo(FILE *f) {
	assert(f != NULL);
	char buf[1024*1024]; // 三代测序应该不会超过1M.
	int cntReads = 0;
	int maxL = -1;
	int minL = 1000000000;
	while((fgets(buf, sizeof(buf), f)) != NULL) {
		++cntReads;
		if((cntReads&3) == 2) {
			int len = strlen(buf) - 1; // 结尾包括一个换行符
			maxL = maxL>len ?maxL :len;
			minL = minL<len ?minL :len;
		}
	}
	cntReads >>= 2;
	printf("序列条数=        %d\n", cntReads);
	printf("序列长度最大值=  %d\n", maxL);
	printf("序列长度最小值=  %d\n", minL);
}

int getReadsInfo_main(int argc, char *argv[]) {
	if(argc == 1) {
		getReadsInfoUsage();
		return 0;
	}
	int c;
	FILE *f = NULL;
	while((c=getopt(argc, argv, "f:")) > -1){
		if(c == 'f') {
			f = fopen(optarg, "r");
			if(f == NULL) {
				fprintf(stderr, "fail to open FASTQ file\n");
				return 1;
			}
		} else {
			fprintf(stderr, "option error!\n");
			getReadsInfoUsage();
			return 1;
		}
	}
	getReadsInfo(f);
	if(f != NULL) fclose(f);
	return 0;
}