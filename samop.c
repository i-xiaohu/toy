//
// Created by 63175 on 2019/7/31.
//

#include <assert.h>
#include <getopt.h>
#include "samop.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

SAMOpt_t samGetOpt(kstring_t *a, SAMOpAux_t *aux) {
	kstring_t *buf = &aux->optBuf; buf->l = 0;
	kstring_t *bs = &aux->bigStr;
	SAMOpt_t res; memset(&res, 0, sizeof(res));
	int i, stop = 0;
	for(i=0; i<=a->l; ++i) {
		if(i==a->l || a->s[i]==':') {
			++stop;
			if(stop == 1) { // TAG
				res.tag = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs); // 注意必须要手动添加'\0'作为结束符。
			} else if(stop == 2) { // TYPE
				res.type = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			} else if(stop == 3) { // VALUE
				res.value = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			}
			buf->l = 0; // 清空l，但内存m不会变。
		} else {
			kputc(a->s[i], buf);
		}
	}
	return res;
}

SAMNode_t samGetNode(char *s, SAMOpAux_t *aux) {
	SAMNode_t a; memset(&a, 0, sizeof(SAMNode_t));
	kstring_t *buf = &aux->nodeBuf; buf->l = 0;
	kstring_t *bs = &aux->bigStr;
	SAMOpt_v *opts = &aux->opts;
	a.optL = opts->n;
	int i, stop=0, len=(int)strlen(s);
	for(i=0; i<=len; ++i) {
		if(i==len || s[i]=='\t') {
			++stop;
			if(stop == 1) { // qname
				a.qname = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			} else if(stop == 2) { // flag
				a.flag = kSGetNumber(buf);
			} else if(stop == 3) { // rname
				a.rname = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			} else if(stop == 4) { // pos
				a.pos = kSGetNumber(buf);
			} else if(stop == 5) { // MapQ
				a.mapq = kSGetNumber(buf);
			} else if(stop == 6) { // CIGAR
				a.cigar = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			} else if(stop == 7) { // ref name of next read
				a.rnext = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			} else if(stop == 8) { // position of next read
				a.pnext = kSGetNumber(buf);
			} else if(stop == 9) { // observed template length
				a.tlen = kSGetNumber(buf);
			} else if(stop == 10) { // Seq
				a.seq = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			} else if(stop == 11) { // Qual
				a.qual = bs->l;
				kputsn(buf->s, buf->l, bs); kputc('\0', bs);
			} else if(stop >= 12) { // Options
				kv_push(SAMOpt_t, *opts, samGetOpt(buf, aux));
			}
			buf->l = 0;
		} else {
			kputc(s[i], buf);
		}
	}
	a.optR = opts->n;
	return a;
}

void samHandleHeader(char *s, SAMHeader_t *h, SAMOpAux_t *aux) {
	kstring_t *buf = &aux->headerBuf; buf->l = 0;
	kstring_t *bs = &aux->bigStr;
	int i, len=(int)strlen(s), type = -1;
	/* 获得type, 一共5种 */
	for(i=0; i<len; ++i) {
		if(s[i]=='\t') {
			if(strcmp(buf->s, "@HD") == 0) {
				type = 1;
			} else if(strcmp(buf->s, "@SQ") == 0) {
				type = 2;
			} else if(strcmp(buf->s, "@RG") == 0) {
				type = 3;
			} else if(strcmp(buf->s, "@PG") == 0) {
				type = 4;
			} else {
				type = 5;
			}
			buf->l = 0;
			++i;
			break;
		} else {
			kputc(s[i], buf);
		}
	}
	assert(type != -1);
	if(type != 5) {
		/* 以TAB为分隔符，末尾一定有个换行符，可能有多个TAG:VALUE */
		SAMHeaderTV_v *tvs;
		if(type == 1) {
			tvs = &h->HD;
		} else if(type == 2) { // SQ和RG都可能有多行
			tvs = kv_pushp(SAMHeaderTV_v, h->SQs);
		} else if(type == 3) {
			tvs = kv_pushp(SAMHeaderTV_v, h->RGs);
		} else {
			tvs = &h->PG;
		}
		kv_init(*tvs); // 这里可以清空，每次处理一行的信息。
		for(; i<=len; ++i) {
			if(i==len || s[i]=='\t') {
				SAMHeaderTV_t tv; memset(&tv, 0, sizeof(tv));
				// 所有的tag都是2个字符
				tv.tag = bs->l;
				kputsn(buf->s, 2, bs);  kputc('\0', bs);
				tv.value = bs->l;
				kputsn(buf->s+3, buf->l-3, bs); kputc('\0', bs);
				kv_push(SAMHeaderTV_t, *tvs, tv);
				buf->l = 0;
			} else {
				kputc(s[i], buf);
			}
		}
	} else {
		/* 一行注释信息 */
		for(; i<len; ++i) {
			kputc(s[i], buf);
		}
		h->CO = bs->l;
		kputsn(buf->s, buf->l, bs);  kputc('\0', bs);
	}
}
void samShowHeader(SAMHeader_t *h, SAMOpAux_t *aux) {
	kstring_t *bs = &aux->bigStr;
	int i, j;
	printf("\n");
	if(h->HD.n != 0) {
		printf("----------@HD----------\n");
		for(i=0; i<h->HD.n; ++i) {
			SAMHeaderTV_t *tv = &h->HD.a[i];
			if(i > 0) printf("\t");
			printf("%s:%s", bs->s+tv->tag, bs->s+tv->value);
		}
		printf("\n");
	}
	if(h->SQs.n != 0) {
		printf("----------@SQ----------\n");
		for(i=0; i<h->SQs.n; ++i) {
			SAMHeaderTV_v *tvs = &h->SQs.a[i];
			for(j=0; j<tvs->n; ++j) {
				SAMHeaderTV_t *tv = &tvs->a[j];
				if(j > 0) printf("\t");
				printf("%s:%s", bs->s+tv->tag, bs->s+tv->value);
			}
			printf("\n");
		}
	}
	if(h->RGs.n != 0) {
		printf("----------@RG----------\n");
		for(i=0; i<h->RGs.n; ++i) {
			SAMHeaderTV_v *tvs = &h->RGs.a[i];
			for(j=0; j<tvs->n; ++j) {
				SAMHeaderTV_t *tv = &tvs->a[j];
				if(j > 0) printf("\t");
				printf("%s:%s", bs->s+tv->tag, bs->s+tv->value);
			}
			printf("\n");
		}
	}
	if(h->PG.n != 0) {
		printf("----------@PG----------\n");
		for(i=0; i<h->PG.n; ++i) {
			SAMHeaderTV_t *tv = &h->PG.a[i];
			if(i > 0) printf("\t");
			printf("%s:%s", bs->s+tv->tag, bs->s+tv->value);
		}
		printf("\n");
	}
	if(h->CO != -1) {
		printf("----------@CO----------\n");
		printf("%s\n", bs->s+h->CO);
	}
	printf("\n");
}

void samInitAux(SAMOpAux_t *aux) {
	memset(&aux->bigStr, 0, sizeof(aux->bigStr));
	memset(&aux->headerBuf, 0, sizeof(aux->headerBuf));
	kv_init(aux->opts);
	memset(&aux->optBuf, 0, sizeof(aux->optBuf));
	memset(&aux->nodeBuf, 0, sizeof(aux->nodeBuf));
}

void samFreeAux(SAMOpAux_t *aux) {
	free(aux->bigStr.s);
	free(aux->headerBuf.s);
	free(aux->opts.a);
	free(aux->optBuf.s);
	free(aux->nodeBuf.s);
}

SAMInfo_t samGetInfo(FILE *f) {
	char buf[1024 * 1024];
	SAMInfo_t info;
	SAMHeader_t *header = &info.header; memset(header, 0, sizeof(*header));
	header->CO = -1;
	SAMNode_v *nodes = &info.rec1; kv_init(*nodes);
	SAMOpAux_t *aux = &info.aux; samInitAux(aux);
	int cntHeader = 0;
	size_t bytes = 0;
	while(fgets(buf, sizeof(buf), f) != NULL) {
		size_t bLen = strlen(buf);
		if(buf[bLen-1] == '\n') { // 去掉最后的换行符，fgets一定会把换行符读进来。
			buf[bLen-1] = '\0';
		}
		bytes += bLen;
		if(buf[0] == '@') {
			++cntHeader;
			samHandleHeader(buf, header, aux);
		} else {
			kv_push(SAMNode_t, *nodes, samGetNode(buf, aux));
			if(nodes->n == 1) {
				samShowHeader(header, aux);
			} else {
				if(nodes->n%1000000 == 0) {
					printf("[%s] has processed %ld SAM nodes, %ld bytes\n", __func__, nodes->n, bytes);
				}
			}
		}
	}
	printf("[%s] get %d SAM header information, %ld bytes\n", __func__, cntHeader, bytes);
	printf("[%s] get %ld SAM nodes\n", __func__, nodes->n);
	return info;
}

void samDestroyInfo(SAMInfo_t *info) {
	// 首部信息没有实际内存，SAM节点内部也没有实际内存，辅助结构中有所有的实际内存

	/* 释放辅助结构 */
	SAMOpAux_t *aux = &info->aux;
	samFreeAux(aux);

	/* 释放头部数组 */
	int i;
	SAMHeader_t *h = &info->header;
	free(h->HD.a);
	for(i=0; i<h->SQs.n; ++i) { // SQ和RG为二维数组
		SAMHeaderTV_v *tvs = &h->SQs.a[i];
		free(tvs->a);
	}
	free(h->SQs.a);
	for(i=0; i<h->RGs.n; ++i) {
		SAMHeaderTV_v *tvs = &h->RGs.a[i];
		free(tvs->a);
	}
	free(h->RGs.a);
	free(h->PG.a);

	/* 释放节点数组 */
	SAMNode_v *v = &info->rec1;
	free(v->a);
}

void samShowSingle(SAMNode_t *q, SAMOpAux_t *aux, FILE *f) {
	if(f == NULL) f = stdout;
	kstring_t *bs = &aux->bigStr;
	fprintf(f, "qname:    %s\n", bs->s+q->qname);
	fprintf(f, "flag:     %d\n", q->flag);
	if((q->flag&samFMulSeg) != 0) {
		fprintf(f, "    template having multiple segments in sequencing\n");
	}
	if((q->flag&samFSegAli) != 0) {
		fprintf(f, "    each segment properly aligned according to the aligner\n");
	}
	if((q->flag&samFUnmap) != 0) {
		fprintf(f, "    segment unmapped\n");
	}
	if((q->flag&samFNxtSeq) != 0) {
		fprintf(f, "    next segment in the template unmapped\n");
	}
	if((q->flag&samFRC) != 0) {
		fprintf(f, "    SEQ being reverse complemented\n");
	}
	if((q->flag&samFNRC) != 0) {
		fprintf(f, "    SEQ of the next segment in the template being reverse complemented\n");
	}
	if((q->flag&samFFst) != 0) {
		fprintf(f, "    the first segment in the template\n");
	}
	if((q->flag&samFLst) != 0) {
		fprintf(f, "   the last segment in the template\n");
	}
	if((q->flag&samFSecAli) != 0) {
		fprintf(f, "    secondary alignment\n");
	}
	if((q->flag&samFNPF) != 0) {
		fprintf(f, "    not passing filters, such as platform/vendor quality controls\n");
	}
	if((q->flag&samFPCR) != 0) {
		fprintf(f, "    PCR or optical duplicate\n");
	}
	if((q->flag&samFSupAli) != 0) {
		fprintf(f, "    supplementary alignment\n");
	}
	fprintf(f, "rname:    %s\n", bs->s+q->rname);
	fprintf(f, "pos:      %d\n", q->pos);
	fprintf(f, "mapq:     %d\n", q->mapq);
	fprintf(f, "CIGAR:    %s\n", bs->s+q->cigar);
	fprintf(f, "rnext:    %s\n", bs->s+q->rnext);
	fprintf(f, "pnext:    %d\n", q->pnext);
	fprintf(f, "tlen:     %d\n", q->tlen);
	fprintf(f, "seq:      %s\n", bs->s+q->seq);
	fprintf(f, "qual:     %s\n", bs->s+q->qual);

	size_t i, optL=q->optL, optR=q->optR;
	SAMOpt_v *opts = &aux->opts;
	fprintf(f, "options:  %ld\n", optR-optL);
	for(i=optL; i<optR; ++i) {
		SAMOpt_t *opt = &opts->a[i];
		fprintf(f, "    %s:%s:%s\n", bs->s+opt->tag, bs->s+opt->type, bs->s+opt->value);
	}
	fprintf(f, "\n");
}

inline void samDumpNode(SAMNode_t *sam, SAMOpAux_t *aux, FILE *fo) {
	kstring_t *bs = &aux->bigStr;
	SAMOpt_v *opts = &aux->opts;
	fprintf(fo, "%s", bs->s+sam->qname);
	fprintf(fo, "\t%d", sam->flag);
	fprintf(fo, "\t%s", bs->s+sam->rname);
	fprintf(fo, "\t%d", sam->pos);
	fprintf(fo, "\t%d", sam->mapq);
	fprintf(fo, "\t%s", bs->s+sam->cigar);
	fprintf(fo, "\t%s", bs->s+sam->rnext);
	fprintf(fo, "\t%d", sam->pnext);
	fprintf(fo, "\t%d", sam->tlen);
	fprintf(fo, "\t%s", bs->s+sam->seq);
	fprintf(fo, "\t%s", bs->s+sam->qual);
	size_t j;
	for(j=sam->optL; j<sam->optR; ++j) {
		SAMOpt_t *opt = &opts->a[j];
		fprintf(fo, "\t%s:%s:%s", bs->s+opt->tag, bs->s+opt->type, bs->s+opt->value);
	}
	fprintf(fo, "\n");
}

void samDump(SAMInfo_t *info, FILE *fo) {
	SAMHeader_t *h = &info->header;
	SAMOpAux_t *aux = &info->aux;
	kstring_t *bs = &aux->bigStr;
	int i, j;
	/* 输出头部信息 */
	if(h->HD.n > 0) {
		fprintf(fo, "@HD");
		for(i=0; i<h->HD.n; ++i) {
			SAMHeaderTV_t *tv = &h->HD.a[i];
			fprintf(fo, "\t%s:%s", bs->s+tv->tag, bs->s+tv->value);
		}
		fprintf(fo, "\n");
	}
	if(h->SQs.n > 0) {
		for(i=0; i<h->SQs.n; ++i) {
			SAMHeaderTV_v *tvs = &h->SQs.a[i];
			fprintf(fo, "@SQ");
			for(j=0; j<tvs->n; ++j) {
				SAMHeaderTV_t *tv = &tvs->a[j];
				fprintf(fo, "\t%s:%s", bs->s+tv->tag, bs->s+tv->value);
			}
			fprintf(fo, "\n");
		}
	}
	if(h->RGs.n > 0) {
		for(i=0; i<h->RGs.n; ++i) {
			SAMHeaderTV_v *tvs = &h->RGs.a[i];
			fprintf(fo, "@RG");
			for(j=0; j<tvs->n; ++j) {
				SAMHeaderTV_t *tv = &tvs->a[j];
				fprintf(fo, "\t%s:%s", bs->s+tv->tag, bs->s+tv->value);
			}
			fprintf(fo, "\n");
		}
	}
	if(h->PG.n > 0) {
		fprintf(fo, "@PG");
		for(i=0; i<h->PG.n; ++i) {
			SAMHeaderTV_t *tv = &h->PG.a[i];
			fprintf(fo, "\t%s:%s", bs->s+tv->tag, bs->s+tv->value);
		}
		fprintf(fo, "\n");
	}
	if(h->CO != -1) {
		fprintf(fo, "@CO");
		fprintf(fo, "\t%s\n", bs->s+h->CO);
	}
	/* 输出SAMNode信息 */
	SAMNode_v *v = &info->rec1;
	for(i=0; i<v->n; ++i) {
		SAMNode_t *sam = &v->a[i];
		samDumpNode(sam, aux, fo);
		fflush(fo); // 其实我认为这步没有必要，因为不是格式化输出占用的缓冲区过大导致的段错误
	}
}

/* 所谓deduplicate就是处理一条序列对应多个SAM records的情况
 * 从目前的数据来看，只需要把带supplumentary alignment和secondary alignment从SAM records中去掉即可
 * 把剩下的SAM records输出即可以作为SAM2SFQ的输入，生成SFQ格式文件，这是BWA-MEZ的输入
 * multi SAMs' reads的被输出到dup文件中去 */
int samDedup(SAMInfo_t *info, FILE *dup) {
	SAMNode_v *v = &info->rec1;
	char *SO = samGetHeaderValue(info, "SO");
	if(SO==NULL || strcmp(SO, "queryname")!=0) {
		fprintf(stderr, "[%s] you have sort SAM file by queryname\n", __func__);
		return 1;
	}
	SAMOpAux_t *aux = &info->aux;
	printf("[%s] raw SAM has %ld nodes\n", __func__, v->n);
	size_t i, j;
	kstring_t *bs = &aux->bigStr;
	/* 开始去冗余的过程，将对info->sams进行修改 */
	size_t oldN = v->n;
	v->n = 0;
	size_t cntMulti = 0;
	for(i=0; i<oldN; ++i) {
		SAMNode_t *si = &v->a[i];
		SAMNode_t *ps = si; // 拥有完整SEQ的比对
		// 找所有跟sam[i]有相同queryname的sam[j]
		for(j=i+1; j<oldN; ++j) {
			SAMNode_t *sj = &v->a[j];
			if(strcmp(bs->s+si->qname, bs->s+sj->qname) == 0) {
				if(dup != NULL) {
					if(j == i+1) {
						++cntMulti;
						/* 第一条比对si一般为主要比对，不带其他flag标识 */
						fprintf(dup, "multiple SAMs' read %ld\n", cntMulti);
						samShowSingle(si, aux, dup);
					}
					/* 剩余的比对sj，一般为非主要比对（前提是排序得当） */
					samShowSingle(sj, aux, dup);
				}
				if((sj->flag&samFSecAli)==0 && (sj->flag&samFSupAli)==0) {
					ps = sj;
				}
			} else {
				break;
			}
		}
		i = j-1;
		kv_push(SAMNode_t, *v,*ps);
	}
	printf("[%s] multiple SAMs' reads number = %ld\n", __func__, cntMulti);
	/* 所谓的去冗余，其实是为了生成sfq文件，真实压缩软件排序不具备这么好的效果 */
	printf("[%s] deduplicated SAM has %ld nodes\n", __func__, v->n);
	return 0;
}

typedef struct {
	size_t all; // 一条序列可能有多个SAM records.
	size_t secAli; // 一般是由于重复引起的multiple alignment.
	size_t supAli; // chimeric alignment.
	size_t secSup; // 即有multiple alignment, 也有chimeric alignment.
	size_t others; // 其他情况
} mulSAMRecords;

void samGetStatistics(SAMInfo_t *info) {
	SAMNode_v *v = &info->rec1;
	SAMOpAux_t *aux = &info->aux;
	int cntUnmapped = 0;
	size_t i, j;
	kstring_t *bs = &aux->bigStr;
	mulSAMRecords mulSR; memset(&mulSR, 0, sizeof(mulSR));
	for(i=0; i<v->n; ++i) {
		SAMNode_t *si = &v->a[i];
		if((si->flag&samFUnmap) != 0) {
			++cntUnmapped; // 未比对上的序列不可能有Multiple SAM records
		}
		int fMulSAMs = 0, fSecAli = 0, fSupAli = 0;
		for(j=i+1; j<v->n; ++j) {
			SAMNode_t *sj = &v->a[j];
			if(strcmp(bs->s+si->qname, bs->s+sj->qname) != 0) {
				break;
			}
			fMulSAMs = 1;
			if((sj->flag&samFSecAli) != 0) fSecAli = 1;
			if((sj->flag&samFSupAli) != 0) fSupAli = 1;
		}
		if(fMulSAMs == 1) {
			++mulSR.all;
			if(fSecAli == 1 && fSupAli == 0) ++mulSR.secAli; // 这条序列为重复序列
			else if(fSecAli == 0 && fSupAli == 1) ++mulSR.supAli; // 这条序列为chimeric read
			else if(fSecAli == 1 && fSupAli == 1) ++mulSR.secSup; // 即重复，又是chimeric, 这种情况我认为应该不会出现
			else ++mulSR.others; // 由其他情况导致了该序列有多条SAM records
		}
		i = j-1;
	}
	printf("[%s] SAM records                     = %ld\n", __func__, v->n);
	printf("[%s] unmapped reads                  = %d\n", __func__, cntUnmapped);
	printf("[%s] reads with multiple SAM records = %ld\n", __func__, mulSR.all);
	printf("[%s]     multiple mapping      = %ld\n", __func__, mulSR.secAli);
	printf("[%s]     chimeric alignments   = %ld\n", __func__, mulSR.supAli);
	printf("[%s]     multiple and chimeric = %ld\n", __func__, mulSR.secSup);
	printf("[%s]     others                = %ld\n", __func__, mulSR.others);
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    samop handle SAM file\n");
	fprintf(stderr, "Usage:      samop [options] <in.sam>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -s INT        show single SAM record\n");
	fprintf(stderr, "            -o STR        dump SAM records into output file\n");
	fprintf(stderr, "            -d            remove chimeric reads，only keep the soft clipped read\n");
	fprintf(stderr, "            -p            remove chimeric reads，only keep the soft clipped read\n");
	fprintf(stderr, "            -D STR        dump chimeric reads into output file\n");
	fprintf(stderr, "\n");
	return 1;
}

int samop_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c, dedup=0, show=-1;
	FILE *f_sam, *f_out, *f_dup;
	f_sam = f_out = f_dup = NULL;
	while((c=getopt(argc, argv, "s:o:dD:")) > -1) {
		if(c == 'o') {
			f_out = fopen(optarg, "w");
			if(f_out == NULL) {
				fprintf(stderr, "fail to open output SAM file %s\n", optarg);
				return 1;
			}
		} else if(c == 's') {
			show = atoi(optarg);
		} else if(c == 'd') {
			dedup = 1;
		} else if(c == 'D') {
			f_dup = fopen(optarg, "w");
			if(f_dup == NULL) {
				fprintf(stderr, "file to open Duplicated file %s\n", optarg);
				return 1;
			}
		} else {
			fprintf(stderr, "option error!\n");
			return usage();
		}
	}
	if(f_dup != NULL) { dedup = 1; }
	if(optind+1 >= argc || optind+2 < argc) {
		return usage();
	}
	f_sam = fopen(argv[optind], "r");
	if(f_sam == NULL) { fprintf(stderr, "[%s] failed to open [%s]\n", __func__, argv[optind]); return 1; }

	SAMInfo_t info = samGetInfo(f_sam);
	samGetStatistics(&info); // 输出unmapped的统计信息
	SAMNode_v *sams = &info.rec1;
	if(show>=0 && show<sams->n) {
		samShowSingle(&sams->a[show], &info.aux, stdout);
	}
	if(dedup == 1) {
		samDedup(&info, f_dup);
	}
	if(f_out != NULL) {
		samDump(&info, f_out);
	}

	/* 释放SAMInfo所占空间 */
	samDestroyInfo(&info);

	/* 关闭所有文件 */
	if(f_sam != NULL) fclose(f_sam);
	if(f_out!= NULL) fclose(f_out);
	if(f_dup != NULL) fclose(f_dup);
	return 0;
}