//
// Created by 63175 on 2019/11/22.
//

#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <assert.h>

#include "kseq.h"
#include "sync_pe.h"
#include "kstring.h"
#include "ksort.h"
#include "progress.h"

KSEQ_DECLARE(gzFile)
void *kopen(const char *fn, int *_fd);
int kclose(void *a);

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    sync-pe   adjust reads2 to fit reads1 QNAME order\n");
	fprintf(stderr, "Usage:      sync-pe   [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -1 STR        reads1 input FASTQ file\n");
	fprintf(stderr, "            -2 STR        reads2 input FASTQ file\n");
	fprintf(stderr, "            -o STR        reodered reads2 output FASTQ file\n");
	fprintf(stderr, "\n");
	return 1;
}

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}

static inline char *dupkstring(const kstring_t *str, int dupempty)
{
	char *s = (str->l > 0 || dupempty)? malloc(str->l + 1) : NULL;
	if (!s) return NULL;

	memcpy(s, str->s, str->l);
	s[str->l] = '\0';
	return s;
}

static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{
	s->name = dupkstring(&ks->name, 1);
	s->comment = dupkstring(&ks->comment, 0);
	s->seq = dupkstring(&ks->seq, 1);
	s->qual = dupkstring(&ks->qual, 0);
	s->l_seq = ks->seq.l;
}

bseq1_t *bseq_read1(int chunk_size, int *n_, long *bytes_, void *ks_){
	kseq_t *ks = (kseq_t*)ks_;
	int size = 0, m = 0, n = 0;
	bseq1_t *seqs = NULL;
	long bytes = 0;
	while ((bytes = kseq_read(ks)) >= 0) {
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		trim_readno(&ks->name);
		kseq2bseq1(ks, &seqs[n]);
		seqs[n].id = n;
		size += seqs[n++].l_seq;
		if (chunk_size != -1 && size >= chunk_size && (n&1) == 0) break;
	}
	*bytes_ = bytes;
	*n_ = n;

	return seqs;
}

typedef struct {
	int id;
	uint64_t qn_hash;
} hseq_t;

static inline uint64_t string_hash(const char *s) {
	uint64_t res = 0;
	int i;
	for(i = 0; s[i]; ++i) {
		res *= 64;
		res += s[i];
	}
	return res;
}

#define hseq_lt(a, b) ((a).qn_hash < (b).qn_hash)
KSORT_INIT(hseq, hseq_t, hseq_lt)

static void do_sync(kseq_t *ks1, kseq_t *ks2, FILE *fo) {
	int i, j, n1, n2; long bytes1 = 0, bytes2 = 0;
	bseq1_t *seqs1 = bseq_read1(-1, &n1, &bytes1, ks1);
	bseq1_t *seqs2 = bseq_read1(-1, &n2, &bytes2, ks2);
	assert(n1 == n2);
	hseq_t *hs2 = malloc(n2 * sizeof(hseq_t));
	for(i = 0; i < n2; ++i) {
		hs2[i].id = i;
		hs2[i].qn_hash = string_hash(seqs2[i].name);
	}
	fprintf(stderr, "sort by qname-hash value for seqs2...\n");
	ks_introsort(hseq, n2, hs2);
	fprintf(stderr, "output adjusted reads2...\n");
	progress_t bar;
	progress_init(&bar, "", 100, PROGRESS_CHR_STYLE);
	for(i = 0; i < n1; ++i) {
		uint64_t h = string_hash(seqs1[i].name);
		int l = 0, r = n2-1, mid, end = -1;
		while(l <= r) {
			mid = (l + r) >> 1;
			if(hs2[mid].qn_hash >= h) {
				end = mid;
				r = mid - 1;
			} else {
				l = mid + 1;
			}
		}
		if(end == -1 || hs2[end].qn_hash != h) {
			fprintf(stderr, "name1 = %s\n", seqs1[i].name);
			if(end != -1) fprintf(stderr, "name2 = %s\n", seqs2[hs2[end].id].name);
		}
		assert(end != -1 && hs2[end].qn_hash == h);
		bseq1_t *res = NULL;
		for(j = end; j < n2; ++j) {
			bseq1_t *p = &seqs2[hs2[j].id];
			if(!strcmp(seqs1[i].name, p->name)) {
				res = p;
				break;
			}
		}
		assert(res != NULL);
		fputs("@", fo); fputs(res->name, fo); fputs("/2\n", fo);
		fputs(res->seq, fo); fputs("\n", fo);
		fputs("+\n", fo);
		fputs(res->qual, fo); fputs("\n", fo);
		progress_show(&bar, 1.0 * (i+1) / n1);
	}
	progress_destroy(&bar);
	free(hs2);
	for(i = 0; i < n1; ++i) {
		bseq1_t *p = &seqs1[i];
		free(p->name); free(p->comment); free(p->seq); free(p->qual);
		p = &seqs2[i];
		free(p->name); free(p->comment); free(p->seq); free(p->qual);
	}
	free(seqs1);
	free(seqs2);
	fprintf(stderr, "\n");
	fprintf(stderr, "done!\n");
}

int sync_pe_main(int argc, char *argv[]) {
	if (argc == 1) {
		return usage();
	}
	int fd1 = 0, fd2 = 0;
	gzFile fp1 = 0, fp2 = 0;
	void *ko1 = NULL, *ko2 = NULL;
	FILE *fo = NULL;

	const char *opts = "1:2:o:";
	int c;
	while((c=getopt(argc, argv, opts)) > -1){
		if(c == '1') {
			ko1 = kopen(optarg, &fd1);
		} else if(c == '2') {
			ko2 = kopen(optarg, &fd2);
		} else if(c == 'o') {
			fo = fopen(optarg, "w");
		} else {
			return usage();
		}
	}
	if(ko1 == NULL || ko2 == NULL || fo == NULL) {
		fprintf(stderr, "[E::%s] fail to open file\n", __func__);
		return 1;
	}
	fp1 = gzdopen(fd1, "r"); fp2 = gzdopen(fd2, "r");
	kseq_t *ks1 = kseq_init(fp1), *ks2 = kseq_init(fp2);
	do_sync(ks1, ks2, fo);
	kseq_destroy(ks1); gzclose(fp1); kclose(ko1);
	kseq_destroy(ks2); gzclose(fp2); kclose(ko2);
	return 0;
}
