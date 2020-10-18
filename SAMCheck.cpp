#include <iostream>
#include <string>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
#include <vector>
#include <unistd.h>
#include <assert.h>
#include "samop.h"

using namespace std;

/*
 * 检查两个SAM文件的一致性
 */

void check2SAM(vector<SAMNode_t> *zv, vector<SAMNode_t> *mv) {
    int i, n=zv->size();
    int z1m0 = 0; // MEZ比对成功，MEM比对失败
    int z1m1 = 0; // MEZ比对成功，MEM比对成功
    int z0m0 = 0; // MEZ比对失败，MEM比对失败
    int z0m1 = 0; // MEZ比对失败，MEM比对成功
    int cntNM = 0; // MEZ与MEM都比对成功时，编辑距离的一致率
    int cntPos = 0; // MEZ与MEM都比对成功时，
    for(i=0; i<n; ++i) {
        SAMNode_t *zs = &zv->at(i);
        SAMNode_t *ms = &mv->at(i);
        assert(zs->qname == ms->qname);
        if(zs->rname!="*" && ms->rname=="*") {
            ++z1m0;
        } else if(zs->rname!="*" && ms->rname!="*") {
            ++z1m1;
            int zNM = samGetOptNM(&zs->options);
            assert(zNM != -1); // 比对上的序列NM不应该等于-1
            int mNM = samGetOptNM(&ms->options);
            assert(zNM != -1);
            if(zNM == mNM) {
                ++cntNM;
            }
            if(zs->pos == ms->pos) {
                ++cntPos;
            }
        } else if(zs->rname=="*" && ms->rname=="*") {
            ++z0m0;
        } else if(zs->rname=="*" && ms->rname!="*") {
            ++z0m1;
        }
    }
    printf("MEZ0    MEM0    %d\n", z0m0);
    printf("MEZ0    MEM1    %d\n", z0m1);
    printf("MEZ1    MEM0    %d\n", z1m0);
    printf("MEZ1    MEM1    %d\n", z1m1);
    printf("NM一致率 = %.3f\n", 100.0*cntNM/z1m1);
    printf("Pos一致率= %.3f\n", 100.0*cntPos/z1m1);
}

void getSAMStatistics(vector<SAMNode_t> *v) {
    int cntRNameStar = 0;
    int cntUnmapped = 0;
    for(int i=0; i<v->size(); ++i) {
        SAMNode_t *sam = &v->at(i);
        if(sam->rname == "*") {
            ++cntRNameStar;
        }
        if((sam->flag&samFUnmap) != 0) {
            ++cntUnmapped;
        }
    }
    printf("RName==*  的数量 = %d\n", cntRNameStar);
    printf("unmapped  的数量 = %d\n", cntUnmapped);
}

void SAMCheckUsage() {
    cerr << endl;
    cerr << "Program:    SAMCheck checks the identical rate of two SAM files" << endl;
    cerr << "Usage:      SAMCheck [options]"  << endl;
    cerr << endl;
    cerr << "            -z STR        MEZ-SAM file name" << endl;
    cerr << "            -m STR        MEM-SAM file name" << endl;
    cerr << endl;
}

int SAMCheck_main(int argc, char *argv[]) {
    if(argc == 1) {
    	SAMCheckUsage();
        return 0;
    }
    int c;
    FILE *fMEZ, *fMEM;
    fMEZ = fMEM = NULL;
    while((c=getopt(argc, argv, "m:z:")) > -1){
        if(c == 'z'){
            fMEZ = fopen(optarg, "r");
            if(fMEZ == NULL){
                cerr << "fail to open MEZ file\n";
                return 1;
            }
        } else if(c == 'm'){
            fMEM = fopen(optarg, "r");
            if(fMEM == NULL){
                cerr << "fail to open MEM file\n";
                return 1;
            }
        } else{
            cerr << "option error!" << endl;
            SAMCheckUsage();
            return 1;
        }
    }
    vector<SAMNode_t> mezSAMs, memSAMs;
    samGetNodes(fMEZ, &mezSAMs);
    // getSAMStatistics(&mezSAMs);

    samGetNodes(fMEM, &memSAMs);
    // getSAMStatistics(&memSAMs);

    if(mezSAMs.size() == memSAMs.size()) {
        check2SAM(&mezSAMs, &memSAMs);
    } else {
        cerr << "两个SAM文件大小不一致，可能需要去冗余操作" << endl;
    }
    fclose(fMEZ);
    fclose(fMEM);
    return 0;
}