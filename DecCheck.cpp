/*************************************************************************
	> File Name: check.cpp
	> Author: ixiaohu
	> Mail:
	> Created Time: 2018年12月02日 星期日 16时00分59秒
 ************************************************************************/

/*
 * 该程序检查解压缩是否正确（检查两个fastq文件的一致性）
 */

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <unistd.h>
#include <stdlib.h>
#include <string>

using namespace std;

int DecCheck_main(int argc, char *argv[]) {
    if(argc == 1){
        fprintf(stderr, "check whether the decompression is right.\n");
        fprintf(stderr, "-q [fastq filename]\n");
        fprintf(stderr, "-a [fasta filename]\n");
        fprintf(stderr, "-d [decompressed filename]\n");
        fprintf(stderr, "-n [reads_n]\n");
        return 0;
    }
    string *a, *b;
    int c, n=0, flag=0;
    FILE *fd, *f2;
    fd = f2 = NULL;
    while((c=getopt(argc, argv, "q:a:d:n:")) >= 0){
        if(c == 'd') {
            fd = fopen(optarg, "r");
        }else if(c == 'a'){
            f2 = fopen(optarg, "r");
        }else if(c == 'q'){
            flag = 1;
            f2 = fopen(optarg, "r");
        }else if(c == 'n'){
            n = atoi(optarg);
        }else {
            fprintf(stderr, "option error!\n");
            return 1;
        }
    }
    if(fd==NULL || f2==NULL){
        fprintf(stderr, "file error!\n");
        return 1;
    }
    a = new string[n+10];
    b = new string[n+10];
    char buf[1024];
    int na=0, nb=0, i;
    for(i=0; i<n; ++i){
        fgets(buf, sizeof(buf), fd);
        a[na++] = buf;
    }
    for(i=0; i<n; ++i){
        fgets(buf, sizeof(buf), f2);
        if(flag)  fgets(buf, sizeof(buf), f2);
        b[nb++] = buf;
        if(flag)  fgets(buf, sizeof(buf), f2);
        if(flag)  fgets(buf, sizeof(buf), f2);
    }
    if(na != nb){
        fprintf(stderr, "two files have different size.\n");
        return 1;
    }
    sort(a, a+na);
    sort(b, b+na);
    int wCnt = 0;
    for(int i=0; i<na; ++i){
        if(a[i] != b[i]){
            fprintf(stderr, "seq[%d], there's one read error!\n", i+1);
            fprintf(stderr, "\t");
            std::cerr << a[i];
            fprintf(stderr, "\t");
            std::cerr << b[i];
            fprintf(stderr, "\n");
            ++wCnt;
            if(wCnt > 10){
                fprintf(stderr, "10 more errors, program is over!\n");
                return 0;
            }
        }
    }
    if(wCnt == 0){
        fprintf(stderr, "no error happens!\n");
    }
    return 0;
}
