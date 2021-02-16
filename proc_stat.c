//
// Created by ixiaohu on 2020/5/7.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#define PNAME_LINE  1
#define VMSIZE_LINE 13
#define VMRSS_LINE  17

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: proc-stat [options] <ID|name> \n");
	fprintf(stderr, "Note: Stats from /proc/pid/status\n");
	fprintf(stderr, "    -p  [INT] Process ID\n");
	fprintf(stderr, "    -n  [STR] Process name\n");
	fprintf(stderr, "    -d  [INT] Delay time [3s]");
	fprintf(stderr, "\n");
	return 1;
}

static char* get_pname(unsigned int pid) {
	char file_name[512], line_buff[512];
	memset(file_name, 0, sizeof(file_name));
	memset(line_buff, 0, sizeof(line_buff));
	FILE *fd;
	sprintf(file_name,"/proc/%d/status",pid);

	fd = fopen(file_name, "r");
	if(fd == NULL) return NULL;
	int i;
	for (i = 0; i < PNAME_LINE; i++){
		fgets(line_buff, sizeof(line_buff), fd);
	}
	fclose(fd);

	char name[64], pname[64];
	sscanf(line_buff,"%s %s", name, pname);
	char *ret = strdup(pname);
	return ret;
}

static long get_VM(unsigned int pid){
	char file_name[512], line_buff[512];
	memset(file_name, 0, sizeof(file_name));
	memset(line_buff, 0, sizeof(line_buff));
	FILE *fd;
	sprintf(file_name,"/proc/%d/status",pid);

	fd = fopen(file_name, "r");
	if(fd == NULL) return -1;
	int i;
	for (i = 0; i < VMSIZE_LINE; i++){
		fgets(line_buff, sizeof(line_buff), fd);
	}
	fclose(fd);

	char name[64];
	long ret;
	sscanf(line_buff,"%s %ld", name, &ret);
	return ret;
}

static long get_RAM(unsigned int pid) {
	char file_name[512], line_buff[512];
	memset(file_name, 0, sizeof(file_name));
	memset(line_buff, 0, sizeof(line_buff));
	FILE *fd;
	sprintf(file_name,"/proc/%d/status",pid);

	fd = fopen(file_name,"r");
	if(fd == NULL) return -1;
	int i;
	for (i = 0; i < VMRSS_LINE; i++) {
		fgets(line_buff, sizeof(line_buff), fd);
	}
	fclose(fd);

	char name[64];
	long ret;
	sscanf(line_buff,"%s %ld", name, &ret);
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


int proc_stat_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c, delay = 3;

	unsigned int pid = 0;
	char *pname = NULL;
	while((c=getopt(argc, argv, "p:n:d:")) >= 0) {
		if (c == 'd') {
			delay = atoi(optarg);
		}
		else if (c == 'p') {
			pid = atoi(optarg);
			pname = get_pname(pid);
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
		now_VM = get_VM(pid);
		if(now_VM == -1) {
			fprintf(stderr, "Process [%d, %s] finished.\n", pid, pname);
			break;
		}
		max_VM = now_VM > max_VM ? now_VM : max_VM;
		now_RAM = get_RAM(pid);
		max_RAM = now_RAM > max_RAM ? now_RAM : max_RAM;
		fprintf(stderr, "RAM: %.2f MB\n", 1.0 * now_RAM / 1024);
		fprintf(stderr, "VM:  %.2f MB\n", 1.0 * now_VM / 1024);
		sleep(delay);
	}
	fprintf(stderr, "Max_RAM: %.2f MB\n", 1.0 * max_RAM / 1024);
	fprintf(stderr, "Max_VM:  %.2f MB\n", 1.0 * max_VM / 1024);

	free(pname);
	return 0;
}
