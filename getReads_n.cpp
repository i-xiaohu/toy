#include <iostream>
#include <cstdio>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <cstring>

#include "utils.h"
using namespace std;

/*
 * 获得fastq | fastq文件中的序列条数
 */

int getReads_n_main(int argc, char *argv[]) {
    if(argc == 1){
        fprintf(stderr, "this program count reads number in fasta or fastq file\n");
        fprintf(stderr, "-q [fastq filename]\n");
        fprintf(stderr, "-a [fasta filename]\n");
        return 0;
    }
    int c, isQ=0;
    long long n = 0, L = 0;
    FILE *f;
    f = NULL;
    while((c=getopt(argc, argv, "a:q:")) > -1) {
        if(c == 'q'){
            f = fopen(optarg, "r");
            isQ = 1;
        }else if(c == 'a'){
            f = fopen(optarg, "r");
        }else {
            fprintf(stderr, "option error!\n");
            return 1;
        }
    }
    if(f == NULL) {
        fprintf(stderr, "file open failed\n");
        return 1;
    }
    char buf[1024];
    double ctime=cputime(), rtime=realtime();
    while(fgets(buf, sizeof(buf), f)){
        ++n;
        if(isQ){
            fgets(buf, sizeof(buf), f);
            L += strlen(buf)-1;
            fgets(buf, sizeof(buf), f);
            fgets(buf, sizeof(buf), f);
        } else {
            L += strlen(buf)-1;
        }
    }
    printf("reads number = %lld, average length = %f\n", n, 1.0*L/n);
    printf("program finished, (%f cpuTime, %f realTime) elapsed\n", cputime()-ctime, realtime()-rtime);
    return 0;
}