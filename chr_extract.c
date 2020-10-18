//
// Created by 63175 on 2020/4/9.
//

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"

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

KSEQ_DECLARE(gzFile)
void *kopen(const char *fn, int *_fd);
int kclose(void *a);

int chr_ext_main(int argc, char **argv) {
	if(argc == 1) {
		usage();
		return 0;
	}
	int c, show = 0, n = 0;
	FILE *out = NULL;
	int fd = 0; gzFile fp = 0; void *ko = NULL;
	while((c = getopt(argc, argv, "i:sn:o:")) >= 0) {
		if(c == 'i') {
			ko = kopen(optarg, &fd);
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
	if(ko == NULL) {
		fprintf(stderr, "open file failed\n");
		return -1;
	}
	if(n > 0 && out == NULL) return -1;

	fp = gzdopen(fd, "r");
	kseq_t *ks = kseq_init(fp);
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
	kseq_destroy(ks); gzclose(fp); kclose(ko);
	return 0;
}
