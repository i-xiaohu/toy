//
// Created by 63175 on 2019/8/6.
//

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <zlib.h>
#include "sync_pe.h"
#include "kseq.h"

static void usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    reads2fa\n");
	fprintf(stderr, "Usage:      generate fake fastq for reads\n");
	fprintf(stderr, "\n");
}

int reads2fa_main(int argc, char **argv) {
	if(argc == 1) {
		usage();
		return 0;
	}
	char buf[1024];
	int c;
	FILE *fastq = NULL;
	FILE *reads = NULL;
	FILE *fa = NULL;
	while((c = getopt(argc, argv, "i:r:a:")) >= 0) {
		if(c == 'i') {
			fastq = fopen(optarg, "r");
		}
		else if(c == 'r') {
			reads = fopen(optarg, "r");
		}
		else if(c == 'a') {
			fa = fopen(optarg, "w");
		}
	}
	if(fastq == NULL || reads == NULL) return -1;
	while(fgets(buf, sizeof(buf), fastq)) { // qname in fq
		fputs(buf, fa); // output qname to fa
		fgets(buf, sizeof(buf), reads); // bases in reads
		fputs(buf, fa); // output bases to fa
		fgets(buf, sizeof(buf), fastq); // bases in fq
		fgets(buf, sizeof(buf), fastq); // + in fq
		fgets(buf, sizeof(buf), fastq); // qual in fq
	}
	fclose(fa);
	fclose(reads);
	fclose(fastq);
	return 0;
}