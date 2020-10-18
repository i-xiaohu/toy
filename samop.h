//
// Created by 63175 on 2019/7/31.
//

#ifndef TOY_SAMOP_H
#define TOY_SAMOP_H

#include <stdio.h>

#include "kvec.h"
#include "kstring.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/* kstring_t对于大量小字符串来说扩展内存速度过慢，可以直接使用char*指向kstring_t中整体的块内存
 * kstring_t就是字符的容器，里面可以放置多个字符串，以'\0'作为分隔 */

/* 存储SAM节点的选项域 */
typedef struct {
	// size_t表示该字符串在bigStr中起始字符的下标，注意一定要使用size_t，内存占用可能超过4G
	size_t tag;
	size_t type; // 我当时手贱把type写成了int型，结果地址>INT_MAX导致段错误
	size_t value;
} SAMOpt_t;
typedef kvec_t(SAMOpt_t) SAMOpt_v;

/* 辅助结构，提供缓冲区及外部内存。 */
typedef struct {
	kstring_t optBuf; // samGetOpt缓冲区
	kstring_t nodeBuf; // samGetNode缓冲区
	kstring_t headerBuf; // samHandleHeader缓冲区
	kstring_t bigStr; // 这个大字符串存储所有字符串。减少内存分配的次数及碎片内存。
	SAMOpt_v opts;
} SAMOpAux_t;

/* 存储SAM的头部信息数据项 */
typedef struct {
	size_t tag; // 注意不可使用char*指针指向bigStr，因为bigStr的内存也在变
	size_t value;
} SAMHeaderTV_t;
typedef kvec_t(SAMHeaderTV_t) SAMHeaderTV_v;
typedef kvec_t(SAMHeaderTV_v) SAMHeaderTVv_v;

/* 存储SAM的头部信息 */
typedef struct {
	SAMHeaderTV_v HD; // 只能出现一行，且一定在首行
	SAMHeaderTVv_v SQs; // 允许多行，代表染色体编号的顺序
	SAMHeaderTVv_v RGs; // 允许多行，允许无序
	SAMHeaderTV_v PG; // 一行
	size_t CO; // 一行, size_t是bigStr中CO起始字符的下标
} SAMHeader_t;

/* 单个SAM节点，描述一行比对信息 */
typedef struct {
	size_t qname; // bigStr中qname起始字符的下标
	int flag;
	size_t rname; // bigStr中rname起始字符的下标
	int pos;
	int mapq;
	size_t cigar; // bigStr中cigar起始字符的下标
	size_t rnext; // bigStr中rnext起始字符的下标
	int pnext;
	int tlen;
	size_t seq;   // bigStr中seq起始字符的下标
	size_t qual;  // bigStr中qual起始字符的下标
	size_t optL, optR; // 该节点选项域在aux->opts中的始末区间
} SAMNode_t;
typedef kvec_t(SAMNode_t) SAMNode_v;

typedef struct {
	SAMHeader_t header;
	SAMNode_v rec1;
	SAMNode_v rec2;
	SAMOpAux_t aux;
} SAMInfo_t;

#define samFMulSeg 0x1   // template having multiple segments in sequencing
                         // 是否为双端测序
#define samFSegAli 0x2   // each segment properly aligned according to the aligner
                         // 双端序列都被比对上了
#define samFUnmap  0x4   // segment unmapped
#define samFNxtSeq 0x8   // next segment in the template unmapped
#define samFRC     0x10  // SEQ being reverse complemented
#define samFNRC    0x20  // SEQ of the next segment in the template being reverse complemented
#define samFFst    0x40  // the first segment in the template
                         // 双端测序中的第一条
#define samFLst    0x80  // the last segment in the template
                         // 双端测序中的第二条
#define samFSecAli 0x100 // secondary alignment
                         // multiple mapping中出现的选项
#define samFNPF    0x200 // not passing filters, such as platform/vendor quality controls
#define samFPCR    0x400 // PCR or optical duplicate
#define samFSupAli 0x800 // supplementary alignment
                         // chimeric alignment中出现的选项

/* 从SAM文件f中读入所有的SAM相关信息 */
SAMInfo_t samGetInfo(FILE *f);

/* 释放所有的SAMNode节点以及头部信息 */
void samDestroyInfo(SAMInfo_t *v);

/* 展示一条比对信息 */
void samShowSingle(SAMNode_t *q, SAMOpAux_t *aux, FILE *f);


/* 对SAM文件进行去冗余，保留一条序列多个比对位置中的一条 */
int samDedup(SAMInfo_t *info, FILE *fo);

/* 将所有的SAM信息再输出到文件里 */
void samDump(SAMInfo_t *info, FILE *fo);

/* 取得编辑距离 */
static inline int samGetOptNM(SAMNode_t *sam, SAMOpAux_t *aux) {
	size_t i, j, l=sam->optL, r=sam->optR; // 注意这里一定是size_t，可能会超INT_MAX
	kstring_t *bs = &aux->bigStr;
	SAMOpt_v *opts = &aux->opts;
	for(i=l; i<r; ++i) {
		SAMOpt_t *opt = &opts->a[i];
		if(strcmp(bs->s+opt->tag, "NM") == 0) {
			int res = 0;
			for(j=opt->value; bs->s[j]!='\0'; ++j) {
				res *= 10;
				res += bs->s[j]-'0';
			}
			return res;
		}
	}
	return -1;
}

/* 获得对应@HD中的相应键值 */
static inline char* samGetHeaderValue(SAMInfo_t *info, const char *tag) {
	SAMHeaderTV_v *v = &info->header.HD;
	kstring_t *bs = &info->aux.bigStr;
	int i;
	for(i=0; i<v->n; ++i) {
		SAMHeaderTV_t *tv = &v->a[i];
		if(strcmp(bs->s+tv->tag, tag) == 0) {
			return bs->s+tv->value;
		}
	}
	return NULL;
}

/* rname==*的数量和umapped的数量 */
void samGetStatistics(SAMInfo_t *info);

#endif //TOY_SAMOP_H
