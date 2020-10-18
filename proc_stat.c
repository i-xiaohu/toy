//
// Created by 63175 on 2020/5/7.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include "utils.h"

#define VMSIZE_LINE 13
#define VMRSS_LINE 17
#define TIME_LMT 600

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:      proc-stat [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            --pid         process ID\n");
	fprintf(stderr, "            --pname       process name\n");
	fprintf(stderr, "\n");
	return 1;
}

static unsigned int get_proc_virtualmem(unsigned int pid){
	char file_name[64]={0};
	FILE *fd;
	char line_buff[512]={0};
	sprintf(file_name,"/proc/%d/status",pid);

	fd = fopen(file_name, "r");
	if(fd == NULL){
		return 0;
	}

	char name[64];
	int i, vmsize;
	for (i = 0; i < VMSIZE_LINE - 1; i++){
		fgets(line_buff, sizeof(line_buff), fd);
	}

	fgets(line_buff, sizeof(line_buff), fd);
	sscanf(line_buff,"%s %d",name,&vmsize);
	fclose(fd);

	return vmsize;
}

static unsigned int get_proc_mem(unsigned int pid) {
	char file_name[64]; memset(file_name, 0, sizeof(file_name));
	FILE *fd;
	char line_buff[512]; memset(line_buff, 0, sizeof(line_buff));
	sprintf(file_name,"/proc/%d/status",pid);

	fd = fopen(file_name,"r");
	if(fd == NULL){
		return 0;
	}

	char name[64];
	int i, vmrss;
	for (i = 0; i < VMRSS_LINE - 1; i++){
		fgets(line_buff, sizeof(line_buff), fd);
	}

	fgets(line_buff,sizeof(line_buff),fd);
	sscanf(line_buff,"%s %d",name, &vmrss);
	fclose(fd);

	return vmrss;
}

static int get_pid(const char* process_name, const char* user) {
	if(user == NULL){
		user = getlogin();
	}

	char cmd[512];
	if (user){
		sprintf(cmd, "pgrep %s -u %s", process_name, user);
	}
	FILE *pstr = popen(cmd,"r");
	if(pstr == NULL){
		return 0;
	}
	char buff[512];
	memset(buff, 0, sizeof(buff));
	if(fgets(buff, 512, pstr) == NULL){
		return 0;
	}
	return atoi(buff);
}


int proc_stat_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c, lo_index;
	const char *short_opts = "";
	const struct option long_opts[] = {
			// {name, has_arg, flag, val}
			{"pid", required_argument, NULL, 0},
			{"pname", required_argument, NULL, 0},
			{NULL, 0, NULL, 0}
	};

	unsigned int pid = 0;
	while((c=getopt_long(argc, argv, short_opts, long_opts, &lo_index)) >= 0) {
		switch (c) {
			case 0:
				if(!strcmp(long_opts[lo_index].name, "pid")) {
					pid = atoi(optarg);
				}
				else if(!strcmp(long_opts[lo_index].name, "pname")) {
					pid = get_pid(optarg, NULL);
				}
				break;
			default : return usage();
		}
	}
	fprintf(stderr, "Process ID: %d\n", pid);
	if(pid > 0) {
		double rtime = realtime();
		int max_mem = 0, now_mem;
		int max_vmem = 0, now_vmem;
		while((now_vmem = get_proc_virtualmem(pid)) > 0) {
			if(realtime() - rtime > TIME_LMT) {
				break;
			}
			max_vmem = now_vmem > max_vmem ?now_vmem :max_vmem;
			now_mem = get_proc_mem(pid);
			max_mem = now_mem > max_mem ?now_mem :max_mem;
		}
		fprintf(stderr, "Max memory(RAM): %.2f M\n", 1.0 * max_mem / 1024);
		fprintf(stderr, "Max virtual memory(RAM): %.2f M\n", 1.0 * max_vmem / 1024);
	}
	return 0;
}
