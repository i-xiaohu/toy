#include <stdio.h>
#include <string.h>

#include "utils.h"

int hfastq_main(int argc, char **argv);
int samop_main(int argc, char *argv[]);
int sam2sortfq_main(int argc, char *argv[]);
int reads2fa_main(int argc, char **argv);
int chr_ext_main(int argc, char **argv);
int proc_stat_main(int argc, char *argv[]);
int weval_main(int argc, char *argv[]);

static void show_commands() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: toy (A toy wrote by ixiaohu)\n");
	fprintf(stderr, "Usage:   toy <command> [options]\n\n");
	fprintf(stderr, "Command:\n");
	fprintf(stderr, "    [1] hfastq        handle FASTQ file\n");
	fprintf(stderr, "    [2] samop         handle SAM file\n");
	fprintf(stderr, "    [3] sam2sfq       SAM -> C-sorted FASTQ file\n");
	fprintf(stderr, "    [4] reads2fa      Add query name to only base reads.\n");
	fprintf(stderr, "    [5] chr-ext       extract one chromosome\n");
	fprintf(stderr, "    [6] proc-stat     status of process\n");
	fprintf(stderr, "    [7] weval         evaluate the SAM file of wgsim reads\n");
	fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
	int i, ret;
	double t_real = realtime();
	if(argc == 1) {
		show_commands();
		return 0;
	}
	if     (strcmp(argv[1], "hfastq") == 0)      ret = hfastq_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "samop") == 0)       ret = samop_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "sam2sfq") == 0)     ret = sam2sortfq_main(argc-1, argv+1);
	else if(strcmp(argv[1], "reads2fa") == 0)    ret = reads2fa_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "chr-ext") == 0)     ret = chr_ext_main(argc-1, argv+1);
	else if(strcmp(argv[1], "proc-stat") == 0)   ret = proc_stat_main(argc-1, argv+1);
	else if(strcmp(argv[1], "weval") == 0)       ret = weval_main(argc-1, argv+1);
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