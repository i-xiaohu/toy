//
// Created by 63175 on 2019/8/6.
//

#include <stdio.h>
#include <unistd.h>
#include <string.h>

void tempUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    temp\n");
	fprintf(stderr, "Usage:      临时写个程序\n");
	fprintf(stderr, "\n");
}

int temp_main(int argc, char *argv[]) {
	if(argc == 1) {
		tempUsage();
		return 0;
	}
	int c;
	FILE *fi=NULL, *fo=NULL;
	while((c=getopt(argc, argv, "i:o:")) > -1){
		if(c == 'i') {
			fi = fopen(optarg, "r");
		} else if(c == 'o') {
			fo = fopen(optarg, "w");
		} else {
			fprintf(stderr, "options error!\n");
			return 1;
		}
	}
	char meta[1024], bases[1024], dis[1024], qual[1024];
	while(1) {
		if(fgets(meta, sizeof(meta), fi) == NULL) {
			break;
		}
		if(fgets(dis, sizeof(dis), fi) == NULL) {
			fprintf(stderr, "[%s] 不够四行\n", __func__);
			return 1;
		}
		if(fgets(bases, sizeof(bases), fi) == NULL) {
			fprintf(stderr, "[%s] 不够四行\n", __func__);
			return 1;
		}
		if(fgets(qual, sizeof(qual), fi) == NULL) {
			fprintf(stderr, "[%s] 不够四行\n", __func__);
			return 1;
		}
		fprintf(fo, "%s", meta);
		fprintf(fo, "%s", bases);
		fprintf(fo, "%s", dis);
		fprintf(fo, "%s", qual);
	}
	if(fi!=NULL) fclose(fi);
	if(fo!=NULL) fclose(fo);
	return 0;
}