#include <stdio.h>
#include <string.h>

#include "utils.h"

int hfastq_main(int argc, char **argv);
int samop_main(int argc, char *argv[]);
int chr_ext_main(int argc, char **argv);
int weval_main(int argc, char *argv[]);
int view_ref_main(int argc, char **argv);
int reorder_qq_main(int argc, char *argv[]);

static void usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: mini_tools (A mini_tools wrote by ixiaohu)\n");
	fprintf(stderr, "Usage:   mini_tools <command> [options]\n\n");
	fprintf(stderr, "Command:\n");
	int cnt = 0;
	fprintf(stderr, "    [%d] hfastq        Handle FASTQ file\n", ++cnt);
	fprintf(stderr, "    [%d] samop         Handle SAM file\n", ++cnt);
	fprintf(stderr, "    [%d] chr-ext       Extract one chromosome\n", ++cnt);
	fprintf(stderr, "    [%d] weval         Evaluate the SAM file of wgsim reads\n", ++cnt);
	fprintf(stderr, "    [%d] view-ref      View reference sequence\n", ++cnt);
	fprintf(stderr, "    [%d] reorder-qq    Reorder query name and quality scores by bases\n", ++cnt);
	fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
	int i, ret;
	double t_real = realtime();
	if(argc == 1) {
		usage();
		return 0;
	}
	if     (strcmp(argv[1], "hfastq") == 0)      ret = hfastq_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "samop") == 0)       ret = samop_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "chr-ext") == 0)     ret = chr_ext_main(argc-1, argv+1);
	else if(strcmp(argv[1], "weval") == 0)       ret = weval_main(argc-1, argv+1);
	else if(strcmp(argv[1], "view-ref") == 0)    ret = view_ref_main(argc-1, argv+1);
	else if(strcmp(argv[1], "reorder-qq") == 0)  ret = reorder_qq_main(argc-1, argv+1);
	else {
		fprintf(stderr, "[%s] unrecognized command '%s'\n", __func__, argv[1]);
		return 1;
	}
	if(ret !=  0) return ret;
	fprintf(stderr, "[%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	return ret;
}