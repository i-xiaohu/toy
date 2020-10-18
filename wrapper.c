//
// Created by 63175 on 2019/10/10.
//

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "utils.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    wrapper executes the shell command --- argv with recording running time.\n");
	fprintf(stderr, "Usage:      wrapper shell command\n");
	fprintf(stderr, "\n");
	return 1;
}

int wrapper_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	double r_time = realtime();
	char cmd[1024]; memset(cmd, 0, sizeof(cmd));
	int i;
	for(i = 1; i < argc; ++i) {
		strcat(cmd, argv[i]);
	}
	fprintf(stderr, "[%s] String execute shell command: %s\n", __func__, cmd);
	system(cmd);
	fprintf(stderr, "[%s] Real time: %.3f, CPU time: %.3f\n", __func__, realtime() - r_time, cputime());
	return 0;
}