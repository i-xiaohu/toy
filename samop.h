//
// Created by 63175 on 2019/1/15.
//

#ifndef SRC_SAMOP_H
#define SRC_SAMOP_H


#include <stdint.h>
#include <string>
#include <vector>

#include "utils.h"
using namespace std;

typedef struct {

} SAMHeader_t;

typedef struct {
	string tag;
	string type;
	string value;
} SAMOpt_t;

typedef struct {
    string qname;
	int flag;
    string rname;
    int pos;
	int mapq;
	string cigar;
	string rnext;
	string pnext;
	int tlen;
    string seq;
    string qual;
    vector<SAMOpt_t> options;
} SAMNode_t;

#define samFMulSeg 0x1   // template having multiple segments in sequencing
#define samFSegAli 0x2   // each segment properly aligned according to the aligner
#define samFUnmap  0x4   // segment unmapped
#define samFNxtSeq 0x8   // next segment in the template unmapped
#define samFRC     0x10  // SEQ being reverse complemented
#define samFNRC    0x20  // SEQ of the next segment in the template being reverse complemented
#define samFFst    0x40  // the first segment in the template
#define samFLst    0x80  // the last segment in the template
#define samFSecAli 0x100 // secondary alignment
#define samFNPF    0x200 // not passing filters, such as platform/vendor quality controls
#define samFPCR    0x400 // PCR or optical duplicate
#define samFSupAli 0x800 // supplementary alignment

/**
 * 从一行比对信息中获得一个结点
 * @param s SAM文件中的一行
 * @param a 存储该行比对信息中的节点
 */
void samGetNode(char *s, SAMNode_t *a);

/**
 * 从SAM文件f中读取所有的SAMNodes存入v中，不返回vector是因为C++ vector的赋值时间与size成线性
 * @param f SAM文件描述符
 * @param v 存储答案的vector指针
 */
void samGetNodes(FILE *f, vector<SAMNode_t> *v);

/**
 * 展示一条比对信息
 * @param q
 */
void samShowSingle(SAMNode_t *q);

/**
 * 对比对信息进行去冗余，仅保留最优比对（编辑距离最小）
 * @param sams
 */
void samDedup(vector<SAMNode_t> *v);

inline int samGetOptNM(vector<SAMOpt_t> *v) {
	for(int i=0; i<v->size(); ++i) {
		SAMOpt_t *opt = &v->at(i);
		if(opt->tag == "NM") {
			return StringToInt(opt->value);
		}
	}
	return -1;
}

#endif //SRC_SAMOP_H
