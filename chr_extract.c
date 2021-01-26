//
// Created by ixiaohu on 2020/4/9.
//

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>

#include "hfastq.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    chr-ext      extract one chromosome\n");
	fprintf(stderr, "Usage:      chr-ext      [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -i STR    input FASTA reference file\n");
	fprintf(stderr, "            -s        show all chromosomes name\n");
	fprintf(stderr, "            -n INT    extracted chromosomes number\n");
	fprintf(stderr, "            -o STR    output FASTA chromosome file\n");
	fprintf(stderr, "\n");
	return 1;
}

int chr_ext_main(int argc, char **argv) {
	if(argc == 1) {
		usage();
		return 0;
	}
	int c, show = 0, n = 0;
	FILE *out = NULL;
	open_fastq_t *of = NULL;
	while((c = getopt(argc, argv, "i:sn:o:")) >= 0) {
		if(c == 'i') {
			of = hfastq_open(optarg);
		}
		else if(c == 's') {
			show = 1;
		}
		else if(c == 'n') {
			n = atoi(optarg);
		}
		else if(c == 'o') {
			out = fopen(optarg, "w");
		}
	}
	if(of == NULL) {
		return -1;
	}
	if(n > 0 && out == NULL) {
		fprintf(stderr, "Please assign a output file.\n");
		return -1;
	}

	kseq_t *ks = of->ks;
	int cnt = 0;
	while(kseq_read(ks) >= 0) {
		cnt++;
		if(show) {
			fprintf(stdout, "No.%d\t%s %s\n", cnt, ks->name.s, ks->comment.s);
		}
		if(cnt == n) {
			fprintf(out, ">%s %s\n", ks->name.s, ks->comment.s);
			fprintf(out, "%s\n", ks->seq.s);
		}
	}
	if(out != NULL) fclose(out);
	hfastq_close(of);
	return 0;
}
