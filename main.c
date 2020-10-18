/* C版本主程序，负责调用其他子程序
 * C
 */
#include <stdio.h>
#include <string.h>

#include "utils.h"

int SAMCheck_main(int argc, char *argv[]);
int getReadsInfo_main(int argc, char *argv[]);
int samop_main(int argc, char *argv[]);
int sam2sfq_main(int argc, char *argv[]);
int sam2pe_sfq_main(int argc, char *argv[]);
int check2Files_main(int argc, char *argv[]);
int temp_main(int argc, char *argv[]);
int multiThread_main(int argc, char *argv[]);
int wrapper_main(int argc, char *argv[]);
int merge_fq_main(int argc, char *argv[]);
int harc2fq_main(int argc, char *argv[]);

void showAllSubCmd() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: toy (A toy wrote by ixiaohu)\n");
	fprintf(stderr, "Usage:   toy <command> [options]\n\n");
	fprintf(stderr, "Command: SAMCheck      check the similarity of two SAM files\n");
	fprintf(stderr, "         getReadsInfo  get reads information of FASTQ file\n");
	fprintf(stderr, "         samop         handle SAM file\n");
	fprintf(stderr, "         sam2sfq       SAM -> sfq file\n");
	fprintf(stderr, "         sam2pe_sfq    SAM -> paired end sfq file\n");
	fprintf(stderr, "         check2Files   check if two files are completely identical\n");
	fprintf(stderr, "         temp          temp program\n");
	fprintf(stderr, "         multiThread   multiple-thread test\n");
	fprintf(stderr, "         wrapper       wrap the specified program with recording running time\n");
	fprintf(stderr, "         merge_fq      merge several FASTQ file\n");
	fprintf(stderr, "         harc2fq       merge HARC decompressed file to FASTQ file\n");
	fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
	int i, ret;
	double t_real = realtime();
	if(argc == 1) {
		showAllSubCmd();
		return 0;
	}
	if(strcmp(argv[1], "SAMCheck") == 0) ret = SAMCheck_main(argc-1, argv+1);
	else if(strcmp(argv[1], "getReadsInfo") == 0) ret = getReadsInfo_main(argc-1, argv+1);
	else if(strcmp(argv[1], "samop") == 0) ret = samop_main(argc-1, argv+1);
	else if(strcmp(argv[1], "sam2sfq") == 0) ret = sam2sfq_main(argc-1, argv+1);
	else if(strcmp(argv[1], "sam2pe-sfq") == 0) ret = sam2pe_sfq_main(argc-1, argv+1);
	else if(strcmp(argv[1], "check2Files") == 0) ret = check2Files_main(argc-1, argv+1);
	else if(strcmp(argv[1], "temp") == 0) ret = temp_main(argc-1, argv+1);
	else if(strcmp(argv[1], "multiThread") == 0) ret = multiThread_main(argc-1, argv+1);
	else if(strcmp(argv[1], "wrapper") == 0) ret = wrapper_main(argc-1, argv+1);
	else if(strcmp(argv[1], "merge_fq") == 0) ret = merge_fq_main(argc-1, argv+1);
	else if(strcmp(argv[1], "harc2fq") == 0) ret = harc2fq_main(argc-1, argv+1);
	else {
		fprintf(stderr, "[%s] unrecognized command '%s'\n", __func__, argv[1]);
		return 1;
	}
	fprintf(stderr, "[%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	return ret;
}