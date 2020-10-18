//
// Created by 63175 on 2019/10/12.
//

#include <stdio.h>
#include <unistd.h>

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    merge_fq (merge several FASTQ files together)\n");
	fprintf(stderr, "Usage:      merge_fq [options] <in1.fastq> [in2.fastq...] \n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -o STR    output merged FASTQ file\n");
	fprintf(stderr, "\n");
	return 1;
}

int merge_fq_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c, i;
	FILE *fo = NULL;
	while((c=getopt(argc, argv, "o:")) > -1) {
		if(c == 'o') {
			fo = fopen(optarg, "w");
		} else {
			return usage();
		}
	}
	fprintf(stderr, "[%s] Starting merging", __func__);
	for(i = optind; i < argc; ++i) {
		fprintf(stderr, " [%s]", argv[i]);
	}
	fprintf(stderr, "\n");

	char buf[1 << 20];
	for(i = optind; i < argc; ++i) {
		FILE *fi = fopen(argv[i], "r");
		if(fi == NULL) {
			fprintf(stderr, "    can't open %s \n", argv[i]);
			return 1;
		}
		while(fgets(buf, sizeof(buf), fi)) {
			fputs(buf, fo);
		}
	}
	fprintf(stderr, "    merging completed\n");
	return 0;
}
