//
// Created by 63175 on 2019/8/5.
//


#include <stdio.h>
#include <unistd.h>
#include <string.h>

void check2FilesUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    check2Files check if \n");
	fprintf(stderr, "Usage:      check2Files [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -1 STR        SAM file1 name\n");
	fprintf(stderr, "            -2 STR        SAM file2 name\n");
	fprintf(stderr, "\n");
}

int check2Files_main(int argc, char *argv[]) {
	if(argc == 1) {
		check2FilesUsage();
		return 0;
	}
	int c;
	FILE *f1 = NULL, *f2 = NULL;
	while((c=getopt(argc, argv, "1:2:")) > -1){
		if(c == '1'){
			f1 = fopen(optarg, "r");
			if(f1 == NULL) {
				fprintf(stderr, "fail to open SAM file1 %s\n", optarg);
				return 1;
			}
		} else if(c == '2') {
			f2 = fopen(optarg, "r");
			if(f2 == NULL) {
				fprintf(stderr, "fail to open SAM file2 %s\n", optarg);
				return 1;
			}
		} else {
			fprintf(stderr, "option error!\n");
			check2FilesUsage();
			return 1;
		}
	}
	char buf1[1024*1024], buf2[1024*1024];
	while(fgets(buf1, sizeof(buf1), f1) != NULL) {
		if(fgets(buf2, sizeof(buf2), f2) != NULL) {
			if(strcmp(buf1, buf2) != 0) {
				fprintf(stderr, "两个文件内容不完全一致\n");
				return 1;
			}
		} else {
			fprintf(stderr, "第二个文件比第一个文件小\n");
		}
	}
	if(fgets(buf2, sizeof(buf2), f2) != NULL) {
		fprintf(stderr, "第一个文件比第二个文件小\n");
		return 1;
	}
	fprintf(stderr, "[%s]: 两个SAM文件完全一致\n", __func__);
	if(f1 != NULL) fclose(f1);
	if(f2 != NULL) fclose(f2);
	return 0;
}