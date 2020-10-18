//
// Created by 63175 on 2019/9/20.
//

#include <stdio.h>

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    sam2pe-sfq\n");
	fprintf(stderr, "Usage:      sam2pe-sfq [options] <in.sam>\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "            -i STR        SAM input file\n");
	fprintf(stderr, "            -1 STR        SFQ file 1\n");
	fprintf(stderr, "            -2 STR        SFQ file 2\n");
	fprintf(stderr, "\n");
	return 1;
}

int sam2pe_sfq_main(int argc, char *argv[]) {
	int ret;
	if(argc == 1) {
		return usage();
	}
	return ret;
}
