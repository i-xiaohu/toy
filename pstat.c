//
// Created by ixiaohu on 2020/5/7.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include "utils.h"
#include "kstring.h"

#define PNAME_LINE  1
#define VMSIZE_LINE 13
#define VMRSS_LINE  16

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: pstat [options] <ID|name> \n");
	fprintf(stderr, "Note: Stats from /proc/pid/status\n");
	fprintf(stderr, "    -p  [INT] Process ID\n");
	fprintf(stderr, "    -n  [STR] Process name\n");
	fprintf(stderr, "    -d  [INT] Delay time [3s]");
	fprintf(stderr, "\n");
	return 1;
}

static char* get_field(int pid, int line_id) {
	char file_name[512], line_buff[512];
	memset(file_name, 0, sizeof(file_name));
	memset(line_buff, 0, sizeof(line_buff));
	FILE *fd;
	sprintf(file_name,"/proc/%d/status",pid);

	fd = fopen(file_name,"r");
	if(fd == NULL) return NULL;
	int i;
	for (i = 0; i < line_id; i++) {
		fgets(line_buff, sizeof(line_buff), fd);
	}
	fclose(fd);

	char name[64], buf[64];
	sscanf(line_buff,"%s %s", name, buf);
	fprintf(stderr, "%s %s\n", name, buf);
	char *ret = strdup(buf);
	return ret;
}

static int get_pid(const char* process_name, const char* user) {
	if(user == NULL) user = getlogin();

	char cmd[512];
	if (user) sprintf(cmd, "pgrep %s -u %s", process_name, user);
	else return -1;
	FILE *pstr = popen(cmd, "r");
	if(pstr == NULL) return -1;

	char buff[512]; memset(buff, 0, sizeof(buff));
	if(fgets(buff, 512, pstr) == NULL) return -1;
	return atoi(buff);
}


int pstat_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	double ctime = cputime(), rtime = realtime();
	kstring_t cmd; memset(&cmd, 0, sizeof(cmd));
	int i;
	for(i = 1; i < argc; i++) {
		kputs(argv[i], &cmd);
		kputc(' ', &cmd);
	}
	system(cmd.s);
	free(cmd.s);
	fprintf(stderr, "%.2f %.2f\n", cputime()-ctime, realtime()-rtime);
	exit(-1);

	int c, delay = 3;

	int pid = 0;
	char *pname = NULL;
	while((c=getopt(argc, argv, "p:n:d:")) >= 0) {
		if (c == 'd') {
			delay = atoi(optarg);
		}
		else if (c == 'p') {
			pid = atoi(optarg);
			pname = get_field(pid, PNAME_LINE);
		}
		else if (c == 'n') {
			pname = strdup(optarg);
			pid = get_pid(pname, "jifahu");
		}
		else return usage();
	}
	if(pid <= 0 || pname == NULL) {
		fprintf(stderr, "Process does not exist.\n");
		return 1;
	}

	fprintf(stderr, "Process status of [%d, %s]\n", pid, pname);
	long max_RAM = 0, now_RAM;
	long max_VM = 0, now_VM;
	while(1) {
		char* fie;
		// Get VM
		fie = get_field(pid, VMSIZE_LINE);
		if(fie == NULL) {
			fprintf(stderr, "Process [%d, %s] finished.\n", pid, pname);
			break;
		}
		now_VM = atol(fie); free(fie);
		max_VM = now_VM > max_VM ? now_VM : max_VM;

		// Get RAM
		fie = get_field(pid, VMRSS_LINE);
		now_RAM = atol(fie); free(fie);
		max_RAM = now_RAM > max_RAM ? now_RAM : max_RAM;
		fprintf(stderr, "RAM: %.2f MB\n", 1.0 * now_RAM / 1024);
		fprintf(stderr, "VM:  %.2f MB\n", 1.0 * now_VM / 1024);
		fprintf(stderr, "\n");
		sleep(delay);
	}
	fprintf(stderr, "Max_RAM: %.2f MB\n", 1.0 * max_RAM / 1024);
	fprintf(stderr, "Max_VM:  %.2f MB\n", 1.0 * max_VM / 1024);

	free(pname);
	return 0;
}
