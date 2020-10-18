//
// Created by 63175 on 2019/7/31.
//
#include <stdio.h>
#include <unistd.h>
#include <assert.h>

#include "samop.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

typedef struct {
	size_t z1m0; // MEZ比对成功，MEM比对失败
	size_t z1m1; // MEZ比对成功，MEM比对成功
	size_t z0m0; // MEZ比对失败，MEM比对失败
	size_t z0m1; // MEZ比对失败，MEM比对成功
	size_t NMZEM; // MEZ/MEM都成功输出一个最优比对，NM不一致的数量
	size_t PosDiff; // MEZ/MEM都成功输出一个最优比对，Pos不一致的数量
	size_t NMZGM; // MEZ-NM > MEM-NM
	size_t NMZLM; // MEZ-NM < MEM-NM
	size_t chiAliZ1M0; // MEZ输出chimeric alignments, 但是MEM没有
	size_t chiAliZ0M1; // MEM输出chimeric alignments, 但是MEZ没有
	size_t chiAliZ1M1; // MEZ与MEM都输出了chimeric alignments, 但是他们的输出情况没有比较
} checkCnt_t;

int check2SAM(SAMInfo_t *zInfo, SAMInfo_t *mInfo, FILE *fo) {
	SAMNode_v *zv = &zInfo->rec1, *mv = &mInfo->rec1;
	SAMOpAux_t *zAux = &zInfo->aux, *mAux = &mInfo->aux;
	kstring_t *zbs = &zAux->bigStr, *mbs = &mAux->bigStr;
	// 对于单个最优mapping，只比较NM
	// 对于chimeric reads，只比较二者是否都为chimeric alignments即可
	char *zSO = samGetHeaderValue(zInfo, "SO");
	if(zSO==NULL || strcmp(zSO, "queryname") != 0) {
		fprintf(stderr, "MEZ-SAM没有按queryname排序\n");
		return 1;
	}
	char *mSO = samGetHeaderValue(mInfo, "SO");
	if(mSO==NULL || strcmp(mSO, "queryname") != 0) {
		fprintf(stderr, "MEM-SAM没有按queryname排序\n");
		return 1;
	}

	checkCnt_t cnt; memset(&cnt, 0, sizeof(cnt));
	size_t mN = mv->n, zN = zv->n;
	size_t mL, mR, zL, zR;
	mL = zL = 0;
	for(; zL<zN && mL<mN; ) {
		SAMNode_t *zsl = &zv->a[zL];
		SAMNode_t *msl = &mv->a[mL];
		assert(strcmp(zbs->s+zsl->qname, mbs->s+msl->qname) == 0);
		for(zR=zL+1; zR < zN; ++zR) {
			SAMNode_t *zsr = &zv->a[zR];
			if(strcmp(zbs->s+zsl->qname, zbs->s+zsr->qname) != 0) {
				break;
			}
		}
		for(mR=mL+1; mR < mN; ++mR) {
			SAMNode_t *msr = &mv->a[mR];
			if(strcmp(mbs->s+msl->qname, mbs->s+msr->qname) != 0) {
				break;
			}
		}
		int zMap = (zsl->flag&samFUnmap)!=0 ?0 :1;
		int mMap = (msl->flag&samFUnmap)!=0 ?0 :1;
		/* 先判断最基本的条件，比对的miss与否 */
		if(zMap==1 && mMap==0) {
			++cnt.z1m0;
		} else if(zMap==0 && mMap==1) {
			++cnt.z0m1;
		} else if(zMap==0 && mMap==0) {
			++cnt.z0m0;
		} else {
			++cnt.z1m1;
			/* MEZ/MEM比对都成功，进行比较评估，需先处理chimeric reads的情况 */
			int zChiAli = zR-zL>1 ?1 :0;
			int mChiAli = mR-mL>1 ?1 :0;
			if(zChiAli==1 && mChiAli==1) {
				++cnt.chiAliZ1M1;
			} else if(zChiAli==1 && mChiAli==0) {
				++cnt.chiAliZ1M0;
			} else if(zChiAli==0 && mChiAli==1) {
				++cnt.chiAliZ0M1;
			} else {
				/* MEZ与MEM都不是chimeric alignments，比较NM即可 */
				int zNM = samGetOptNM(zsl, zAux);
				assert(zNM != -1);
				int mNM = samGetOptNM(msl, mAux);
				assert(mNM != -1);
				if(zNM < mNM) {
					++cnt.NMZLM;
				} else if(zNM > mNM) {
					++cnt.NMZGM;
				} else {
					++cnt.NMZEM;
				}
				if(zsl->pos != msl->pos) {
					++cnt.PosDiff;
				}
			}
		}
		zL = zR;
		mL = mR;
	}
	assert(zL==zN && mL==mN);
	printf("MEZ-unmapped  ~  MEM-unmapped    %ld\n", cnt.z0m0);
	printf("MEZ-unmapped  ~  MEM-mapped      %ld\n", cnt.z0m1);
	printf("MEZ-mapped    ~  MEM-unmapped    %ld\n", cnt.z1m0);
	printf("MEZ-mapped    ~  MEM-mapped      %ld\n", cnt.z1m1);
	printf("    MEZ-chiAli      ~  MEM-chiAli        %ld\n", cnt.chiAliZ1M1);
	printf("    MEZ-chiAli      ~  MEM-not-chiAli    %ld\n", cnt.chiAliZ1M0);
	printf("    MEZ-not-chiAli  ~  MEM-chiAli        %ld\n", cnt.chiAliZ0M1);
	printf("    MEZ-not-chiAli  ~  MEM-not-chiAli\n");
	printf("        MEZ-NM = MEM-NM    %ld\n", cnt.NMZEM);
	printf("        MEZ-NM < MEM-NM    %ld\n", cnt.NMZLM);
	printf("        MEZ-NM > MEM-NM    %ld\n", cnt.NMZGM);
	printf("        MEZ-Pos != MEM-Pos %ld\n", cnt.PosDiff);
	return 0;
}

void SAMCheckUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    SAMCheck checks the identical rate of two SAM files\n");
	fprintf(stderr, "Usage:      SAMCheck [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -z STR        MEZ-SAM file name\n");
	fprintf(stderr, "            -m STR        MEM-SAM file name\n");
	fprintf(stderr, "            -o STR        output NM info into this file\n");
	fprintf(stderr, "\n");
}

int SAMCheck_main(int argc, char *argv[]) {
	if(argc == 1) {
		SAMCheckUsage();
		return 0;
	}
	int c;
	FILE *fMEZ, *fMEM, *fo;
	fMEZ = fMEM = fo = NULL;
	while((c=getopt(argc, argv, "m:z:o:")) > -1){
		if(c == 'z'){
			fMEZ = fopen(optarg, "r");
			if(fMEZ == NULL) {
				fprintf(stderr, "fail to open MEZ file\n");
				return 1;
			}
		} else if(c == 'm'){
			fMEM = fopen(optarg, "r");
			if(fMEM == NULL) {
				fprintf(stderr, "fail to open MEM file\n");
				return 1;
			}
		} else if(c == 'o') {
			fo = fopen(optarg, "w");
			if(fo == NULL) {
				fprintf(stderr, "fail to open output file\n");
				return 1;
			}
		} else {
			fprintf(stderr, "option error!\n");
			SAMCheckUsage();
			return 1;
		}
	}
	SAMInfo_t mezSAMInfo, memSAMInfo;
	mezSAMInfo = samGetInfo(fMEZ);
	// getSAMStatistics(&mezSAMs);

	memSAMInfo = samGetInfo(fMEM);
	// getSAMStatistics(&memSAMs);

	int ret = check2SAM(&mezSAMInfo, &memSAMInfo, fo);
	samDestroyInfo(&mezSAMInfo);
	samDestroyInfo(&memSAMInfo);
	if(fMEZ != NULL) fclose(fMEZ);
	if(fMEM != NULL) fclose(fMEM);
	if(fo != NULL) fclose(fo);
	return ret;
}