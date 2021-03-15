//
// Created by ixiaohu on 2020/5/7.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

#include "utils.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: pstat 'prog args'.\n");
	fprintf(stderr, "Note:  pstat gives the time cost and real memory usage(RAM) of program.\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int i, childpid;
	double rtime = realtime();
	if ((childpid = fork()) == 0) {
		char *cmd_args[argc + 1];
		for(i = 1; i < argc; i++) cmd_args[i-1] = argv[i];
		cmd_args[argc-1] = NULL;
		if (execvp(cmd_args[0], cmd_args) < 0) {
			perror("error on execvp");
			exit(-1);
		}
	} else {
		struct rusage ru;
		int ret, status;
		ret = wait(&status);
		getrusage(RUSAGE_CHILDREN, &ru);
		fprintf(stderr, "\tCommand: \"");
		for(i = 1; i < argc-1; i++) fprintf(stderr, "%s ", argv[i]);
		fprintf(stderr, "%s\"\n", argv[argc-1]);
		fprintf(stderr, "\tMAX_rss: %ld kbytes\n", ru.ru_maxrss);
		double cpu_sec = ru.ru_utime.tv_sec + ru.ru_stime.tv_sec + 1e-6 * (ru.ru_utime.tv_usec + ru.ru_stime.tv_usec);
		fprintf(stderr, "\tCPU_seconds: %.3f\n",cpu_sec );
		fprintf(stderr, "\tElapsed_seconds: %.3f\n", realtime() - rtime);
		fprintf(stderr, "\tParent_exit: child_pid=%d, wait_return=%d, child status=%d\n", childpid, ret, status);
	}
	return 0;
}
