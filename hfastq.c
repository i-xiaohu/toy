//
// Created by ixiaohu on 2019/7/31.
//

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <getopt.h>
#include <zlib.h>

#include "kseq.h"
#include "utils.h"
#include "progress.h"
#include "table.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    hfastq       get the FASTQ file information.\n");
	fprintf(stderr, "Usage:      hfastq       [options] <data1.fq> [data2.fq]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -K   INT       chunk size of inputting once [100M]\n");
	fprintf(stderr, "\n");
	return 1;
}

KSEQ_DECLARE(gzFile)
void *kopen(const char *fn, int *_fd);
int kclose(void *a);

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

/** Reading FASTQ sequences from bwa.h */
typedef struct {
	int l_seq, id;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

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

/** sequences statistics */
typedef struct {
	int seqs_n;
	int min_l, max_l, reads_N;
	size_t bases_N, bases_all;
	int single_max_N, single_min_N;
} seqs_info_t;

typedef struct {
	int chunk_size;
	int show_bar;
} opt_t;

static opt_t* opt_init() {
	opt_t *a = calloc(1, sizeof(opt_t));
	a->chunk_size = 100*1000*1000;
	a->show_bar = 1;
	return a;
}

typedef struct {
	opt_t *opt;
	kseq_t *ks1, *ks2;
	long f1_size, f2_size, has_read_bytes;
	long n_processed;
	progress_t bar;
	seqs_info_t seqs1_info, seqs2_info;
} ktp_aux_t;

static void aux_init(ktp_aux_t *a) {
	memset(a, 0, sizeof(*a));
	progress_init(&a->bar, "", 100, PROGRESS_CHR_STYLE);
	a->seqs1_info.min_l = INT_MAX;
	a->seqs2_info.min_l = INT_MAX;
	a->seqs1_info.single_min_N = INT_MAX;
	a->seqs2_info.single_min_N = INT_MAX;
}

static void aux_detroy(ktp_aux_t *a) {
	progress_destroy(&a->bar);
}

typedef struct {
	int n_seqs1, n_seqs2;
	bseq1_t *seqs1, *seqs2;
} ktp_data_t;

static int compare_pe_qname(const char *s1, const char *s2) {
	int i;
	for(i = 0; s1[i] != '/' && s1[i] != '\0'; ++i) {
		if(s1[i] != s2[i]) { return 0; }
	}
	return 1;
}

static void seqs_info(seqs_info_t *s, int n, bseq1_t *seqs) {
	s->seqs_n += n;
	int i, j;
	for(i = 0; i < n; i++) {
		const bseq1_t *p = &seqs[i];
		s->min_l = s->min_l < p->l_seq ?s->min_l :p->l_seq;
		s->max_l = s->max_l > p->l_seq ?s->max_l :p->l_seq;
		s->bases_all += p->l_seq;
		int N = 0;
		for(j = 0; j < p->l_seq; j++) {
			N += (p->seq[j] == 'N');
		}
		if(N > 0) {
			s->reads_N++;
			s->bases_N += N;
			s->single_min_N = s->single_min_N < N ?s->single_min_N :N;
			s->single_max_N = s->single_max_N > N ?s->single_max_N :N;
		}
	}
}

static void *process(void *shared, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	opt_t *opt = aux->opt;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if(step == 0) {
		// Input reads
		ktp_data_t *ret;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs1 = bseq_read1(opt->chunk_size, &ret->n_seqs1, &aux->has_read_bytes, aux->ks1);
		if(aux->ks2) {
			ret->seqs2 = bseq_read1(opt->chunk_size, &ret->n_seqs2, &aux->has_read_bytes, aux->ks2);
			assert(ret->n_seqs1 == ret->n_seqs2);
		} else {
			ret->seqs2 = NULL;
			ret->n_seqs2 = 0;
		}
		if(ret->seqs1 == NULL) {
			free(ret);
			return NULL;
		}
		return ret;
	} else if(step == 1) {
		// Process reads
		aux->n_processed += data->n_seqs1;
		seqs_info(&aux->seqs1_info, data->n_seqs1, data->seqs1);
		if(aux->ks2) {
			seqs_info(&aux->seqs2_info, data->n_seqs2, data->seqs2);
		}
		aux->has_read_bytes = aux->has_read_bytes <= 0 ?(aux->f1_size + aux->f2_size) :aux->has_read_bytes;

		if(opt->show_bar) {
			progress_show(&aux->bar, 1.0 * aux->has_read_bytes / (aux->f1_size + aux->f2_size));
			if(aux->has_read_bytes == aux->f1_size + aux->f2_size) {
				fprintf(stderr, "\n");
			}
		} else {
			fprintf(stderr, "[%s] Read %d sequences.\n", __func__, data->n_seqs1 + data->n_seqs2);
		}
		return data;
	} else if(step == 2) {
		// Free allocated memory.
		for(i = 0; i < data->n_seqs1; ++i) {
			const bseq1_t *p = &data->seqs1[i];
			free(p->name); free(p->comment); free(p->seq); free(p->qual);
		}
		free(data->seqs1);
		if(aux->ks2) {
			for(i = 0; i < data->n_seqs2; ++i) {
				const bseq1_t *p = &data->seqs2[i];
				free(p->name); free(p->comment); free(p->seq); free(p->qual);
			}
			free(data->seqs2);
		}
		free(data);
		return NULL;
	}
	return NULL;
}

static void output_seqs_info(seqs_info_t *info) {
	info->single_min_N = info->single_min_N == INT_MAX ?0 :info->single_min_N;
	int i, n = 2, m = 6;
	char *h1[25] = {"seqs", "length", "reads N", "bases", "N", "N range"};
	tab_t *t = tab_init(n, m, "seqs statistics", NULL);
	for(i = 1; i <= m; ++i) tab_fill_cell1(t, 1, i, "%s", h1[i-1]);
	tab_fill_cell1(t, 2, 1, "%d", info->seqs_n);
	tab_fill_cell1(t, 2, 2, "%d - %d", info->min_l, info->max_l);
	tab_fill_cell1(t, 2, 3, "%.3f %%", 100.0 * info->reads_N / info->seqs_n);
	tab_fill_cell1(t, 2, 4, "%ld", info->bases_all);
	tab_fill_cell1(t, 2, 5, "%.3f %%", 100.0 * info->bases_N / info->bases_all);
	tab_fill_cell1(t, 2, 6, "%d - %d", info->single_min_N, info->single_max_N);
	tab_display(t); tab_destroy(t);
}

int hfastq_main(int argc, char **argv) {
	if(argc == 1) {
		return usage();
	}
	int c;
	opt_t *opt = opt_init();
	ktp_aux_t aux; aux_init(&aux);
	int fd1 = 0; gzFile fp1 = 0; void *ko1 = NULL;
	int fd2 = 0; gzFile fp2 = 0; void *ko2 = NULL;
	while((c=getopt(argc, argv, "K:")) >= 0) {
		if (c == 'K') {
			opt->chunk_size = atoi(optarg);
		} else {
			return usage();
		}
	}
	if(optind != argc-1 && optind != argc-2) {
		return usage();
	}
	ko1 = kopen(argv[optind], &fd1);
	fp1 = gzdopen(fd1, "r"); aux.ks1 = kseq_init(fp1); aux.f1_size = get_file_size(fd1);
	if(strstr(argv[optind], ".gz")) opt->show_bar = 0;
	if(optind == argc-2) {
		ko2 = kopen(argv[optind + 1], &fd2);
		fp2 = gzdopen(fd2, "r"); aux.ks2 = kseq_init(fp2); aux.f2_size = get_file_size(fd2);
		if(strstr(argv[optind + 1], ".gz")) opt->show_bar = 0;
	}
	aux.opt = opt;

	kt_pipeline(2, process, &aux, 3);

	output_seqs_info(&aux.seqs1_info);
	if(aux.ks2) output_seqs_info(&aux.seqs1_info);

	aux_detroy(&aux);
	kseq_destroy(aux.ks1); gzclose(fp1); kclose(ko1);
	if(ko2 != NULL) {
		kseq_destroy(aux.ks2); gzclose(fp2); kclose(ko2);
	}
	free(opt);
	return 0;
}