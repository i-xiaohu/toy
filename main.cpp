//
// Created by 63175 on 2019/7/28.
//

/*
 * 主程序，负责调用其他子程序
 */

#include <stdio.h>
#include <string.h>

#include "utils.h"
using namespace std;

int calcCov_main(int argc, char *argv[]);
int ctgals_main(int argc, char *argv[]);
int DecCheck_main(int argc, char *argv[]);
int fetch_main(int argc, char *argv[]);
int getReads_n_main(int argc, char *argv[]);
int metaCheck_main(int argc, char *argv[]);
int SAMCheck_main(int argc, char *argv[]);
int splitPE_main(int argc, char *argv[]);

void showAllSubCmd() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: toy (A toy wrote by ixiaohu)\n");
	fprintf(stderr, "Usage:   toy <command> [options]\n\n");
	fprintf(stderr, "Command: calcCov       calculate sequence coverage depth.\n");
	fprintf(stderr, "         ctgals        contig analysis\n");
	fprintf(stderr, "         DecCheck      check if the decompressed file is right\n");
	fprintf(stderr, "         fetch         get several reads(like head)\n");
	fprintf(stderr, "         getReads_n    get reads number\n");
	fprintf(stderr, "         metaCheck     check if the mata is right\n");
	fprintf(stderr, "         SAMCheck      check the similarity of two SAM files\n");
	fprintf(stderr, "         splitPE       split smart paired end file into two files\n");
	fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
	int i, ret;
	double t_real = realtime();
	if(argc == 1) {
		showAllSubCmd();
		return 0;
	}
	if(strcmp(argv[1], "calcCov") == 0) ret = calcCov_main(argc-1, argv+1);
	else if(strcmp(argv[1], "ctgals") == 0) ret = ctgals_main(argc-1, argv+1);
	else if(strcmp(argv[1], "DecCheck") == 0) ret = DecCheck_main(argc-1, argv+1);
	else if(strcmp(argv[1], "fetch") == 0) ret = fetch_main(argc-1, argv+1);
	else if(strcmp(argv[1], "getReads_n") == 0) ret = getReads_n_main(argc-1, argv+1);
	else if(strcmp(argv[1], "metaCheck") == 0) ret = metaCheck_main(argc-1, argv+1);
	else if(strcmp(argv[1], "SAMCheck") == 0) ret = SAMCheck_main(argc-1, argv+1);
	else if(strcmp(argv[1], "splitPE") == 0) ret = splitPE_main(argc-1, argv+1);
	else {
		fprintf(stderr, "[%s] unrecognized command '%s'\n", __func__, argv[1]);
		return 1;
	}
	fprintf(stderr, "[%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
}