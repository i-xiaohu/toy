//
// Created by ixiaohu on 2020/5/7.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

/*
 * execvp + getrusage
 * [bwamem-main] Real time: 180.227 sec; CPU: 3442.331 sec
 *     Max_RAM: 7617864 kbytes
 *     Max_VIRT:  8165064 kbytes
 *     MAX_rss: 7617864 kbytes
 *     CPU Seconds: 3442.351
 * [toy-main] Real time: 180.253 sec; CPU: 83.077 sec
 *
 * /usr/bin/time -v
 * [main] Real time: 219.194 sec; CPU: 4155.078 sec
 *     User time (seconds): 4128.21
 *     System time (seconds): 26.88
 *     Elapsed (wall clock) time (h:mm:ss or m:ss): 3:39.21
 *     Maximum resident set size (kbytes): 30474992
 *
 * Summarize
 *     RSS: Resident set size is the amount of memory the residents(presents) in (real) RAM.
 *     VIRT >= real RAM + swapped memory.
 *     MRSS: Maximum RSS is peak RAM memory usage.
 *     系统调用getrusage()的"ru_maxrss"跟toy pstat(100us每次)访问/proc/pid/status文件中的"VmRSS"最大值一致.
 *
 *     time输出的MRSS大约是其他两个输出的4倍, 据网络资料所言, 这很可能是程序bug.
 *
 *     time和toy pstat的CPU时间和时钟时间都与程序消耗几乎一致.
 *
 * RSS:
 *     RSS是程序实际使用的物理内存
 *     有趣的是, linux上所有进程的RSS加并不等于这些程序实际的资源使用,
 *     RSS多次记录了shared libraries, 所以大于等于实际资源效果.
 */

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: pstat 'cmd args'.\n");
	fprintf(stderr, "    pstat gives the time cost and real memory usage(RAM) of cmd.\n");
	fprintf(stderr, "\n");
	return 1;
}

static long get_field(int pid, const char *filed_name) {
	char file_name[512], line_buff[512];
	memset(file_name, 0, sizeof(file_name));
	memset(line_buff, 0, sizeof(line_buff));
	FILE *fd;
	sprintf(file_name,"/proc/%d/status",pid);

	fd = fopen(file_name,"r");
	if(fd == NULL) return -1;
	long ret = -1;
	while(fgets(line_buff, sizeof(line_buff), fd) != NULL) {
		char id[64], buf[64];
		sscanf(line_buff,"%s %s", id, buf);
		if(!strcmp(id, filed_name)) {
			ret = atol(buf);
			break;
		}
	}
	fclose(fd);
	return ret;
}

int pstat_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int childpid;
	if ((childpid = fork()) == 0) {
		char *cmd_args[argc + 1];
		int i;
		for(i = 1; i < argc; i++) cmd_args[i-1] = argv[i];
		cmd_args[argc-1] = NULL;
		if (execvp(cmd_args[0], cmd_args) < 0) {
			perror("error on execvp");
			exit(-1);
		}
	} else {
		long max_RAM = 0, now_RAM;
		long max_VM = 0, now_VM;
		int pid = getpid();
		while(1) {
			if(get_field(childpid, "PPid:") != pid) {
				fprintf(stderr, "\tChild process [%d] finished.\n", childpid);
				break;
			}

			now_VM = get_field(childpid, "VmSize:");
			if(now_VM == -1) {
				fprintf(stderr, "\tChild process [%d] finished.\n", childpid);
				break;
			}
			max_VM = now_VM > max_VM ? now_VM : max_VM;

			now_RAM = get_field(childpid, "VmRSS:");
			max_RAM = now_RAM > max_RAM ? now_RAM : max_RAM;
			usleep(100);
		}
		struct rusage xxx;

		fprintf(stderr, "\tMax_RAM: %ld kbytes\n", max_RAM);
		fprintf(stderr, "\tMax_VIRT:  %ld kbytes\n", max_VM);
		int ret, status;
		ret = wait(&status);
		getrusage(RUSAGE_CHILDREN, &xxx);
		fprintf(stderr, "\tMAX_rss: %ld kbytes\n", xxx.ru_maxrss);
		fprintf(stderr, "\tCPU Seconds: %.3f\n", 1.0 * (xxx.ru_utime.tv_sec + xxx.ru_stime.tv_sec) + 1e-6 * (xxx.ru_utime.tv_usec+xxx.ru_stime.tv_usec));
		fprintf(stderr, "\tParent exit: child_pid=%d, wait_return=%d, child status=%d\n", childpid, ret, status);
	}
	return 0;
}
