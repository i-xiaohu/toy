//
// Created by 63175 on 2019/10/26.
//

#include <stdio.h>
#include <unistd.h>
#include <string.h>

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    harc2fq merge HARC decompressed file to FASTQ file\n");
	fprintf(stderr, "Usage:      harc2fq [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -i STR        HARC files prefix: [prefix.dna.d], [prefix.id], [prefix.quality]\n");
	fprintf(stderr, "            -o STR        output merged FASTQ file\n");
	fprintf(stderr, "\n");
	return 1;
}

int harc2fq_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c;
	FILE *f_id, *f_qual, *f_seq;
	FILE *f_fq = NULL;
	char *prefix = NULL;
	while((c = getopt(argc, argv, "i:o:")) > -1) {
		if(c == 'i') {
			prefix = optarg;
		} else if(c == 'o') {
			f_fq = fopen(optarg, "w");
		} else {
			return usage();
		}
	}
	int prefix_l = (int)strlen(prefix);
	char fn_id[prefix_l + strlen(".id") + 1]; strcpy(fn_id, prefix); strcat(fn_id, ".id");
	char fn_qual[prefix_l + strlen(".quality") + 1]; strcpy(fn_qual, prefix); strcat(fn_qual, ".quality");
	char fn_seq[prefix_l + strlen(".dna.d") + 1]; strcpy(fn_seq, prefix); strcat(fn_seq, ".dna.d");
	f_id = fopen(fn_id, "r");
	f_qual = fopen(fn_qual, "r");
	f_seq = fopen(fn_seq, "r");
	if(f_id == NULL || f_qual == NULL || f_seq == NULL) {
		fprintf(stderr, "[E] open files failed\n");
		return 1;
	}

	fprintf(stderr, "[%s] merge [%s], [%s], [%s] into FASTQ file\n", __func__, fn_id, fn_seq, fn_qual);

	char buf[300];
	while(fgets(buf, sizeof(buf), f_id) != NULL) {
		fputs(buf, f_fq); // line 1
		fgets(buf, sizeof(buf), f_seq);
		fputs(buf, f_fq); // line 2
		fputs("+\n", f_fq); // line 3
		fgets(buf, sizeof(buf), f_qual);
		fputs(buf, f_fq); // line 4
	}
	return 0;
}
