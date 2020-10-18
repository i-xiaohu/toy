//
// Created by 63175 on 2019/9/1.
//

/* 测试场景：一个线程负责读入，n个线程去处理，不同的数据处理速度不同。
 * 要求：1、读入缓冲区有限，读入线程填满缓冲区后须放弃核心。
 *       2、应满足CPU time ~= real time * n，且不能相差过多，也就是有非常好的并行效果。*/

#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include "kvec.h"

void multiThreadUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    multiple-thread test\n");
	fprintf(stderr, "Usage:      multiThread [options]\n");
	fprintf(stderr, "            -i input FASTQ file\n");
	fprintf(stderr, "            -m per thread memory\n");
	fprintf(stderr, "            -n threads number\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
}

/* 数据单元 */
typedef struct {
	char *meta;
	char *bases;
	char *qual;
} dataNode_t;

/* 每个线程要处理的数据块 */
typedef kvec_t(dataNode_t) dataBlock_t;

/* 读入线程辅助结构 */
struct mtReadAux_t;

/* 读入线程工作结构 */
typedef struct {
	struct mtReadAux_t *aux;
} mtReadWorker_t;

/* 处理线程工作结构 */
typedef struct {
	struct mtReadAux_t *aux;
	int id; // 线程号，用来找对应的数据块号
} mtProWorker_t;

/* 读入线程辅助结构 */
typedef struct mtReadAux_t {
	FILE *fi; // 读入文件
	int processThread_n;  // 处理数据的线程数
	pthread_t rwTId; // 读入线程id
	mtReadWorker_t rw; // 读入线程工作结构
	pthread_t *pwTId; // 处理线程id
	mtProWorker_t *pw; // 处理线程
	size_t perThreadMemory; // 每个线程要处理数据的大小
	size_t lineBufSize; // 行缓冲区大小
	char *lineBuf; // 行缓冲区
	uint8_t *isNull; // 该数据块是否为空
	dataBlock_t *blocks; // 分给线程的数据块
	int eof; // 是否读到文件结尾
	pthread_mutex_t mutex; // 给临界资源visit, eof上锁
	pthread_cond_t cond; // 等待条件：1、读入线程等待visit有空置位置。2、处理线程等待visit不为空
} mtReadAux_t;

inline void mtReadAuxInit(mtReadAux_t *r) {
	r->fi = NULL;
	r->processThread_n = 8;
	r->pwTId = malloc(r->processThread_n * sizeof(pthread_t));
	r->pw = malloc(r->processThread_n * sizeof(mtProWorker_t));
	r->perThreadMemory = (10<<20);
	r->lineBufSize = (10<<10);
	r->lineBuf = malloc(r->lineBufSize * sizeof(char));
	/* 所有isNull置为非空 */
	r->isNull = malloc(r->processThread_n * sizeof(uint8_t));
	memset(r->isNull, 0x1, r->processThread_n * sizeof(uint8_t));
	/* 初始时没有线程占据数据块 */
	/* 初始blocks的{n, m, a}都为空 */
	r->blocks = malloc(r->processThread_n * sizeof(dataBlock_t));
	memset(r->blocks, 0, r->processThread_n * sizeof(dataBlock_t));
	r->eof = 0;
	pthread_mutex_init(&r->mutex, NULL);
	pthread_cond_init(&r->cond, NULL);
}

inline void mtReadAuxDestroy(mtReadAux_t *r) {
	fclose(r->fi);
	free(r->pwTId);
	free(r->pw);
	free(r->lineBuf);
	free(r->isNull);
	free(r->blocks);
	pthread_mutex_destroy(&r->mutex);
	pthread_cond_destroy(&r->cond);
}

/* 读入线程工作函数 */
void *mtReadWorker(void *data) {
	mtReadWorker_t *w = (mtReadWorker_t*)data;
	mtReadAux_t *p = w->aux;

	FILE *fi = p->fi;
	int n = p->processThread_n;
	size_t m = p->perThreadMemory;
	size_t lineBufSize = p->lineBufSize;
	char *lineBuf = p->lineBuf;
	uint8_t *isNull = p->isNull;
	dataBlock_t *blocks = p->blocks;

	for(;;) {
		if(p->eof != 0) { // 只有自己才修改这个变量
			break;
		}
		// 访问临界资源
		pthread_mutex_lock(&p->mutex);
		int i, flag = -1;
		for(;;) {
			for(i=0; i<n; ++i) {
				if(isNull[i] == 1) { // isNull为1，表示此数据块空闲
					flag = i;
					break;
				}
			}
			if(flag != -1) break;
			pthread_cond_wait(&p->cond, &p->mutex);
		}
		pthread_mutex_unlock(&p->mutex);
		assert(flag != -1);
		/* 为空闲位置读入数据 */
		dataBlock_t *block = &blocks[flag];
		block->n = 0; // 数据块buffer都是反复利用的缓冲区
		int lines = 0;
		size_t bytes = 0;
		dataNode_t node; memset(&node, 0, sizeof(node));
		int eof = 0;
		while (1) {
			fgets(lineBuf, lineBufSize, fi);
			eof = feof(fi);
			if(eof != 0) break;
			++lines;
			size_t slen = strlen(lineBuf);
			lineBuf[slen - 1] = '\0';
			--slen;
			bytes += slen + 1;
			if ((lines & 3) == 1) {
				node.meta = malloc(slen + 1);
				strcpy(node.meta, lineBuf);
			} else if ((lines & 3) == 3) {
				node.bases = malloc(slen + 1);
				strcpy(node.bases, lineBuf);
			} else if ((lines & 3) == 0) {
				node.qual = malloc(slen + 1);
				strcpy(node.qual, lineBuf);
				kv_push(dataNode_t, *block, node);
				/* 当读入数据量到达单线程数据块上限时，退出 */
				if (bytes >= m) {
					break;
				}
			}
		}
		fprintf(stderr, "[%s] read %ld sequences into block(%d)\n", __func__, block->n, flag);
		/* 读入完成后，修改isNull及eof，并通知处理线程 */
		pthread_mutex_lock(&p->mutex);
		isNull[flag] = 0;
		p->eof = eof;
		pthread_cond_broadcast(&p->cond);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(NULL);
}

/* 处理线程工作函数 */
void *mtProWorker(void *data) {
	mtProWorker_t *w = (mtProWorker_t*)data;
	mtReadAux_t *p = w->aux;
	int id = w->id;
	uint8_t *isNull = p->isNull;
	for(;;) {
		/* 查看自己的数据块是否有数据 */
		pthread_mutex_lock(&p->mutex);
		int eof = 0;
		for(;;) {
			if(isNull[id]) { // 线程数据块为空
				if(p->eof != 0) { // 到达文件末尾
					eof = 1;
					break;
				} else { // 等待读入线程填满该数据块
					pthread_cond_wait(&p->cond, &p->mutex);
				}
			} else { // 线程数据块不为空
				break;
			}
		}
		pthread_mutex_unlock(&p->mutex);
		if(eof == 1) {
			break;
		}
		/* 将对应数据块取出来处理 */
		dataBlock_t *block = &p->blocks[id];
		int k, n = (rand()%20) * 10000;
		size_t i;
		int ans = 0;
		for(k=0; k<n; ++k) {
			for(i=0; i<block->n; ++i) {
				dataNode_t *a = &block->a[i];
				ans += a->meta[0];
				ans += a->bases[1];
				ans += a->qual[2];
				ans &= 63;
			}
		}
		for(i=0; i<block->n; ++i) {
			dataNode_t *a = &block->a[i];
			free(a->meta);
			free(a->bases);
			free(a->qual);
		}
		fprintf(stderr, "[%s] thread(%d) processed %ld sequences, ans = %d\n", __func__, id, block->n, ans);
		/* 通知其他线程当前数据块已被处理完成 */
		pthread_mutex_lock(&p->mutex);
		isNull[id] = 1;
		pthread_cond_broadcast(&p->cond);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(NULL);
}

int multiThread_main(int argc, char *argv[]) {
	if(argc == 1) {
		multiThreadUsage();
		return 0;
	}
	srand((unsigned int)time(NULL));
	int c;
	mtReadAux_t aux;
	mtReadAuxInit(&aux);
	while((c=getopt(argc, argv, "i:m:n:")) > -1) {
		if(c == 'i') {
			aux.fi = fopen(optarg, "r");
		} else if(c == 'm'){
			char *q;
			aux.perThreadMemory = strtol(optarg, &q, 0);
			if (*q == 'k' || *q == 'K') aux.perThreadMemory <<= 10;
			else if (*q == 'm' || *q == 'M') aux.perThreadMemory <<= 20;
			else if (*q == 'g' || *q == 'G') aux.perThreadMemory <<= 30;
		} else if(c == 'n') {
			aux.processThread_n = atoi(optarg);
			aux.pwTId = realloc(aux.pwTId, aux.processThread_n * sizeof(pthread_t));
			aux.pw = realloc(aux.pw, aux.processThread_n * sizeof(mtProWorker_t));
			aux.isNull = realloc(aux.isNull, aux.processThread_n * sizeof(uint8_t));
			memset(aux.isNull, 0x1, aux.processThread_n * sizeof(uint8_t));
			aux.blocks = realloc(aux.blocks, aux.processThread_n * sizeof(dataBlock_t));
			memset(aux.blocks, 0, aux.processThread_n * sizeof(dataBlock_t));
		} else {
			fprintf(stderr, "option error!\n");
			multiThreadUsage();
			return 1;
		}
	}
	if(aux.fi == NULL) {
		fprintf(stderr, "fail to open FASTQ file\n");
		return 1;
	}
	/* 创建读线程 */
	mtReadWorker_t *rWorker = &aux.rw;
	rWorker->aux = &aux;
	pthread_create(&aux.rwTId, NULL, mtReadWorker, rWorker);

	/* 创建处理线程 */
	int i;
	for(i=0; i<aux.processThread_n; ++i) {
		pthread_t *tid = &aux.pwTId[i];
		mtProWorker_t *pw = &aux.pw[i];
		pw->aux = &aux;
		pw->id = i;
		pthread_create(tid, NULL, mtProWorker, pw);
	}
	pthread_join(aux.rwTId, NULL);
	for(i=0; i<aux.processThread_n; ++i) {
		pthread_join(aux.pwTId[i], NULL);
	}
	mtReadAuxDestroy(&aux);
	return 0;
}
