//
// Created by ixiaohu on 2020/5/7.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>

#include "utils.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: pstat [n|notify] 'prog args'.\n");
	fprintf(stderr, "Note:  pstat gives the time cost and real memory usage(RAM) of program.\n");
	fprintf(stderr, "Parameter notify enables the notification of program status.\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int i, childpid, skip = 1, notify = 0;
	if (!strcmp(argv[1], "n") || !strcmp(argv[1], "notify")) {
		skip = 2;
		notify = 1;
	}

	double rtime = realtime();
	double ctime = 0;
	if ((childpid = fork()) == 0) {
		char *cmd_args[argc - skip + 1];
		for(i = skip; i < argc; i++) cmd_args[i-skip] = argv[i];
		cmd_args[argc-skip] = NULL;
		if (execvp(cmd_args[0], cmd_args) < 0) {
			perror("error on execvp");
			exit(-1);
		}
	} else {
		struct rusage ru;
		int ret, status;
		ret = wait(&status);
		rtime = realtime() - rtime;
		getrusage(RUSAGE_CHILDREN, &ru);
		ctime = ru.ru_utime.tv_sec + ru.ru_stime.tv_sec + 1e-6 * (ru.ru_utime.tv_usec + ru.ru_stime.tv_usec);

		FILE *f = stderr;
		if (notify) {
			char *home_path = getenv("HOME");
			char pstat_path[strlen(home_path) + 16];
			sprintf(pstat_path, "%s/pstat", home_path);
			if (access(pstat_path, 0) == -1) mkdir(pstat_path, 0755);

			char pstat_file[strlen(pstat_path) + 16];
			for (i = 1; i < 1024; i++) {
				sprintf(pstat_file, "%s/%04d.pstat", pstat_path, i);
				if (access(pstat_file, 0) == -1) {
					fprintf(stderr, "\tLog written to notification file %s\n", pstat_file);
					f = fopen(pstat_file, "w");
					break;
				}
			}
			if (f) {
				fprintf(f, "\tCommand: \"");
				for(i = skip; i < argc-1; i++) fprintf(f, "%s ", argv[i]);
				fprintf(f, "%s\"\n", argv[argc-1]);
				fprintf(f, "\tMAX_rss:     ");
				if (ru.ru_maxrss < 1024) fprintf(f, "%ld KB\n", ru.ru_maxrss);
				else if (ru.ru_maxrss < 1024*1024) fprintf(f, "%ld MB\n", ru.ru_maxrss/1024);
				else fprintf(f, "%.1f GB\n", 1.0*ru.ru_maxrss/1024/1024);
				int h = (int)ctime / 3600;
				int m = ((int)ctime - h * 3600) / 60;
				int s = (int)ctime - h * 3600 - m * 60;
				fprintf(f, "\tCPU_time:    %02d:%02d:%02d  i.e.  %.1f sec\n", h, m, s, ctime);
				h = (int)rtime / 3600;
				m = ((int)rtime - h * 3600) / 60;
				s = (int)rtime - h * 3600 - m * 60;
				fprintf(f, "\tWall_clock:  %02d:%02d:%02d  i.e.  %.1f sec\n", h, m, s, rtime);
				fprintf(f, "\tJob_CPU:     %.2f\n", ctime / rtime);
				fprintf(f, "\tParent_exit: child_pid=%d, wait_return=%d, child status=%d\n", childpid, ret, status);
				fclose(f);
			}
		}

		f = stderr;
		fprintf(f, "\tCommand: \"");
		for(i = skip; i < argc-1; i++) fprintf(f, "%s ", argv[i]);
		fprintf(f, "%s\"\n", argv[argc-1]);
		fprintf(f, "\tMAX_rss:     ");
		if (ru.ru_maxrss < 1024) fprintf(f, "%ld KB\n", ru.ru_maxrss);
		else if (ru.ru_maxrss < 1024*1024) fprintf(f, "%ld MB\n", ru.ru_maxrss/1024);
		else fprintf(f, "%.1f GB\n", 1.0*ru.ru_maxrss/1024/1024);
		int h = (int)ctime / 3600;
		int m = ((int)ctime - h * 3600) / 60;
		int s = (int)ctime - h * 3600 - m * 60;
		fprintf(f, "\tCPU_time:    %02d:%02d:%02d  i.e.  %.1f sec\n", h, m, s, ctime);
		h = (int)rtime / 3600;
		m = ((int)rtime - h * 3600) / 60;
		s = (int)rtime - h * 3600 - m * 60;
		fprintf(f, "\tWall_clock:  %02d:%02d:%02d  i.e.  %.1f sec\n", h, m, s, rtime);
		fprintf(f, "\tJob_CPU:     %.2f\n", ctime / rtime);
		fprintf(f, "\tParent_exit: child_pid=%d, wait_return=%d, child status=%d\n", childpid, ret, status);
	}
	return 0;
}
