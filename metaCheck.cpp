/*************************************************************************
	> File Name: metaCheck.cpp
	> Author: jifahu
	> Mail: 
	> Created Time: 2018年12月05日 星期三 10时54分36秒
 ************************************************************************/

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <string>
#include <unistd.h>
using namespace std;
#define N 1000000

/*
 * 检查元数据是否正确
 */

struct MetaNode{
    string meta;
    string bases;
};

inline bool MNCmp(MetaNode A, MetaNode B){
    return A.meta < B.meta;
}

int metaCheck_main(int argc, char *argv[]) {
    if(argc == 1){
        fprintf(stderr, "this program check whether metadata is correct,\n");
        fprintf(stderr, "usage: metaCheck -q [fastq] -m [meta.out] -a [fasta]\n");
        return 0;
    }
    MetaNode a[N], b[N];
    int c;
    FILE *fq, *fm, *fa;
    fq = fm = fa = NULL;
    while((c=getopt(argc, argv, "q:m:a:")) >= 0){
        if(c == 'q'){
            fq = fopen(optarg, "r");
        }else if(c == 'm'){
            fm = fopen(optarg, "r");
        }else if(c == 'a'){
            fa = fopen(optarg, "r");
        }else {
            fprintf(stderr, "usage: metaCheck -q [fastq] -m [meta.out] -a [fasta]\n");
            return 1;
        }
    }
    int na = 0, nb = 0;
    char buf[1024];
    while(fgets(buf, sizeof(buf), fq)){
        a[na].meta = buf;
        fgets(buf, sizeof(buf), fq);
        a[na++].bases = buf;
        fgets(buf, sizeof(buf), fq);
        fgets(buf, sizeof(buf), fq);
    }

    while(fgets(buf, sizeof(buf), fm)){
        b[nb].meta = buf;
        fgets(buf, sizeof(buf), fa);
        b[nb++].bases = buf;
    }
    if(na != nb){
        fprintf(stderr, "error...  na = %d, while nb = %d\n", na, nb);
        return 1;
    }
    fprintf(stderr, "the number of reads is %d\n", na);
    sort(a, a+na, MNCmp);
    sort(b, b+nb, MNCmp);
    for(int i=0; i<na; ++i){
        if(a[i].meta != b[i].meta){
            fprintf(stderr, "after sort, there are two distinct metadata, this is wrong!\n");
            return 1;
        }
        if(a[i].bases != b[i].bases){
            fprintf(stderr, "metadata are same, but not for bases, that means the metadata order is wrong!\n");
            fprintf(stderr, "the %d read in fq is\t", i);
            cout << a[i].bases << endl;
            fprintf(stderr, "            in fa is\t");
            cout << b[i].bases << endl;
            return 1;
        }
    }
    fprintf(stderr, "0 error...\n");
    return 0;
}

