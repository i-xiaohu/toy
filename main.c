#include <stdio.h>
#include <string.h>

#include "utils.h"
#include "table.h"

int hfastq_main(int argc, char **argv);
int samop_main(int argc, char *argv[]);
int check2Files_main(int argc, char *argv[]);
int reads2fa_main(int argc, char **argv);
int multiThread_main(int argc, char *argv[]);
int sam2sortfq_main(int argc, char *argv[]);
int eva2sam_main(int argc, char *argv[]);
int sync_pe_main(int argc, char *argv[]);
int hsfq_main(int argc, char *argv[]);
int spring_reorder_main(int argc, char *argv[]);
int chr_ext_main(int argc, char **argv);
int proc_stat_main(int argc, char *argv[]);

static void show_commands() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: toy (A toy wrote by ixiaohu)\n");
	fprintf(stderr, "Usage:   toy <command> [options]\n\n");
	fprintf(stderr, "Command: eva2sam       evaluate the similarity of two alignment files\n");
	fprintf(stderr, "         hfastq        handle FASTQ file\n");
	fprintf(stderr, "         samop         handle SAM file\n");
	fprintf(stderr, "         sam2sfq       SAM -> C-sorted FASTQ file\n");
	fprintf(stderr, "         check2Files   check if two files are completely identical\n");
	fprintf(stderr, "         reads2fa      \n");
	fprintf(stderr, "         multiThread   multiple-thread test\n");
	fprintf(stderr, "         sync-pe       adjust reads2 to fit reads1 QNAME order\n");
	fprintf(stderr, "         hsfq          handle reordered FASTQ file\n");
	fprintf(stderr, "         reorder       recover reordered FASTQ file by Spring/HARC\n");
	fprintf(stderr, "         chr-ext       extract one chromosome\n");
	fprintf(stderr, "         proc-stat     status of process\n");
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
	else if(strcmp(argv[1], "eva2sam") == 0)     ret = eva2sam_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "check2Files") == 0) ret = check2Files_main(argc-1, argv+1);
	else if(strcmp(argv[1], "reads2fa") == 0)    ret = reads2fa_main(argc - 1, argv + 1);
	else if(strcmp(argv[1], "multiThread") == 0) ret = multiThread_main(argc-1, argv+1);
	else if(strcmp(argv[1], "sam2sfq") == 0)     ret = sam2sortfq_main(argc-1, argv+1);
	else if(strcmp(argv[1], "sync-pe") == 0)     ret = sync_pe_main(argc-1, argv+1);
	else if(strcmp(argv[1], "hsfq") == 0)        ret = hsfq_main(argc-1, argv+1);
	else if(strcmp(argv[1], "reorder") == 0)     ret = spring_reorder_main(argc-1, argv+1);
	else if(strcmp(argv[1], "chr-ext") == 0)     ret = chr_ext_main(argc-1, argv+1);
	else if(strcmp(argv[1], "proc-stat") == 0)   ret = proc_stat_main(argc-1, argv+1);
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