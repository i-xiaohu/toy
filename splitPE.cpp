#include <iostream>
#include <unistd.h>
#include <cstring>
#include <stdio.h>

using namespace std;

/*
 * 把smart paired fastq（智能双端测序文件）分成两个fastq文件。
 */

int splitPE_main(int argc, char *argv[]) {
    if(argc == 1){
        fprintf(stderr, "this program split fastq PE into two file: fn-1.fq and fn-2.fq\n");
        fprintf(stderr, "-i [PE-fastq filename]\n");
        fprintf(stderr, "-p [new filename prefix]\n");
        return 0;
    }
    int c, n1=0, n2=0, _p=0;
    FILE *f = NULL;
    char prefix[128], fn[128];
    while((c=getopt(argc, argv, "i:p:")) > -1) {
        if(c == 'i'){
            strcpy(fn, optarg);
            f = fopen(optarg, "r");
        } else if(c == 'p'){
            _p = 1;
            strcpy(prefix, optarg);
        } else {
            fprintf(stderr, "option error!\n");
            return 1;
        }
    }
    if(!_p) {
        strcpy(prefix, fn);
        size_t lp = strlen(prefix);
        for(size_t i=lp-1; i>=0; --i){
            if(prefix[i] == '.'){
                prefix[i] = '\0';
            }
        }
    }
    if(f == NULL){
        fprintf(stderr, "file error!\n");
        return 1;
    }
    char buf[1024];
    // @SRR554369.1.1 1 length=100
    // @SRR554369.1.2 1 length=100
    char fn1[128], fn2[128];
    strcpy(fn1, prefix); strcat(fn1, "-1.fq");
    strcpy(fn2, prefix); strcat(fn2, "-2.fq");
    FILE *f1, *f2;
    f1 = f2 = NULL;
    f1 = fopen(fn1, "w");
    f2 = fopen(fn2, "w");
    if(f1 == NULL || f2 == NULL){
        fprintf(stderr, "file error!\n");
        return 1;
    }
    while(fgets(buf, sizeof(buf), f)) {
        if(feof(f) || !strcmp(buf,"\n")) break;
        int len = (int)strlen(buf);
        int E = 0;
        for(int i=0; i<len; ++i){
            if(buf[i] == ' '){
                E = buf[i-1]-'0';
                break;
            }
        }
        if(E == 1) {
            fprintf(f1, "%s", buf);
            fgets(buf, sizeof(buf), f); fprintf(f1, "%s", buf);
            fgets(buf, sizeof(buf), f); fprintf(f1, "%s", buf);
            fgets(buf, sizeof(buf), f); fprintf(f1, "%s", buf);
            ++n1;
        } else if(E == 2) {
            fprintf(f2, "%s", buf);
            fgets(buf, sizeof(buf), f); fprintf(f2, "%s", buf);
            fgets(buf, sizeof(buf), f); fprintf(f2, "%s", buf);
            fgets(buf, sizeof(buf), f); fprintf(f2, "%s", buf);
            ++n2;
        } else {
            break;
        }
    }
    fprintf(stderr, "n1 = %d, n2 = %d\n", n1, n2);
    if(n1 != n2){
        fprintf(stderr, "PE-1 != PE-2\n");
        return 1;
    }
    return 0;
}