//
// Created by 63175 on 2019/1/15.
//

/*
 * 该程序计算SAM文件中序列的覆盖深度（平均距离）
 */

#include <iostream>
#include <string>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <unistd.h>
#include "samop.h"

using namespace std;

int calcCov_main(int argc, char *argv[]) {
    if (argc == 1) {
        fprintf(stderr, "this program output the average distance of reads with best AS after sort.\n");
        fprintf(stderr, "-i [SAM filename]\n");
        return 0;
    }
    vector<SAMNode_t> sams;
    int c;
    FILE *fi;
    fi = NULL;
    while((c=getopt(argc, argv, "i:")) > -1) {
        if(c == 'i') {
            fi = fopen(optarg, "r");
            if(fi == NULL){
                cerr << "fail to open SAM file\n";
                return 1;
            }
        }else{
            cerr << "option error!" << endl;
            return 1;
        }
    }
    samGetNodes(fi, &sams);
    int dis=0, close=0;
    for(int i=0;i<sams.size(); ++i) {
        if(i > 0 && sams[i].rname==sams[i-1].rname) {
            dis += sams[i].pos - sams[i-1].pos;
            ++close;
        }
    }
    cout << "the average distance: " << 1.0*dis/close << endl;
    return 0;
}