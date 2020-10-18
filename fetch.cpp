#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

int fetch_main(int argc, char *argv[]){
    if(argc == 1){
        fprintf(stderr, "this program fetch [fastq or fasta reads] from [fq or fa] file\n");
        fprintf(stderr, "-a [fasta filename]\n");
        fprintf(stderr, "-q [fastq filename]\n");
        fprintf(stderr, "-n [reads number]\n");
        return 0;
    }
    int c, isA=1, n=0;
    FILE *f;
    while((c=getopt(argc, argv, "a:q:n:")) >= 0){
        if(c == 'q'){
            isA = 0;
            f = fopen(optarg, "r");
        } else if(c == 'a'){
            f = fopen(optarg, "r");
        } else if(c == 'n'){
            n = atoi(optarg);
        } else {
            fprintf(stderr, "option error!\n");
            return 1;
        }
    }
    char buf[1024];
    int i;
    if(isA){
        for(i=0; i<n; ++i){
            fgets(buf, sizeof(buf), f);
            printf("%s", buf);
        }
    } else {
        n = n<<2;
        for(i=0; i<n; ++i){
            fgets(buf, sizeof(buf), f);
            printf("%s", buf);
        }
    }
    return 0;
}