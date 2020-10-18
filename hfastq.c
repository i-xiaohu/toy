//
// Created by 63175 on 2019/7/31.
//

/*
 * 获得序列相关的信息
 */

#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <getopt.h>
#include <zlib.h>

#include "kseq.h"
#include "utils.h"
#include "sync_pe.h"
#include "progress.h"
#include "table.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/*
 * toy hfastq
 * 1. Offers statistics of FASTQ/FASTA sequences(single and paired end)
 * 2. Simply supports quality control of raw reads
 * 3. Truncates and output partial consecutive reads
 * */

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    hfastq       handle FASTQ file\n");
	fprintf(stderr, "Usage:      hfastq       [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            --i1 STR       input FASTQ file name\n");
	fprintf(stderr, "            --o1 STR       output FASTQ file name\n");
	fprintf(stderr, "            --i2 STR       input FASTQ file name\n");
	fprintf(stderr, "            --o2 STR       output FASTQ file name\n");
	fprintf(stderr, "            --qc INT       perform quality control of raw reads\n");
	fprintf(stderr, "            --io1          not multiple IO thread\n");
	fprintf(stderr, "            -K   INT       chunk size of each thread [10000000]\n");
	fprintf(stderr, "            -t   INT       chunk size of each thread [8]\n");
	fprintf(stderr, "            -l   INT       down sample lower bound\n");
	fprintf(stderr, "            -r   INT       down sample upper bound\n");
	fprintf(stderr, "\n");
	return 1;
}

KSEQ_DECLARE(gzFile)
void *kopen(const char *fn, int *_fd);
int kclose(void *a);

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

typedef struct {
	int seqs_n;
	int min_l, max_l, reads_N;
	size_t bases_N, bases_all;
	int single_max_N, single_min_N;
} seqs_info_t;

typedef struct {
	FILE *input_f1, *input_f2;
	FILE *output_f1, *output_f2;
	int chunk_size;
	int n_threads;
	int qual_ctl;
	int flt_N; // discard reads of N >= INT
} opt_t;

static opt_t* opt_init() {
	opt_t *a = calloc(1, sizeof(opt_t));
	a->chunk_size = 10000000;
	a->n_threads = 8;
	a->flt_N = 10;
	return a;
}

typedef struct {
	opt_t *opt;
	kseq_t *ks1, *ks2;
	long f1_size, f2_size, has_read_bytes;
	long n_processed;
	progress_t bar;
	seqs_info_t seqs1_info, seqs2_info;
	int flt_n; // # read been filter out
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

/************************
 * sequences statistics *
 ************************/
typedef struct {
	seqs_info_t *infos;
	bseq1_t *seqs;
} info_worker_t;

static void seqs_info_worker(void *data, long seq_id, int t_id) {
	info_worker_t *w = (info_worker_t*)data;
	bseq1_t *p = &w->seqs[seq_id];
	seqs_info_t *s = &w->infos[t_id];
	s->min_l = s->min_l < p->l_seq ?s->min_l :p->l_seq;
	s->max_l = s->max_l > p->l_seq ?s->max_l :p->l_seq;
	int i, N = 0;
	s->bases_all += p->l_seq;
	for(i = 0; i < p->l_seq; ++i) {
		N += (p->seq[i] == 'N');
	}
	if(N) {
		++s->reads_N;
		s->bases_N += N;
		s->single_min_N = s->single_min_N < N ?s->single_min_N :N;
		s->single_max_N = s->single_max_N > N ?s->single_max_N :N;
	}
}

static void seqs_info(seqs_info_t *s, int n_threads, int n, bseq1_t *seqs) {
	s->seqs_n += n;
	int i;
	info_worker_t w;
	w.infos = calloc(n_threads, sizeof(seqs_info_t));
	w.seqs = seqs;
	for(i = 0; i < n_threads; ++i) {
		w.infos[i].min_l = INT_MAX;
		w.infos[i].single_min_N = INT_MAX;
	}
	kt_for(n_threads, seqs_info_worker, &w, n);
	for(i = 0; i < n_threads; ++i) {
		const seqs_info_t *p = &w.infos[i];
		s->min_l = s->min_l < p->min_l ?s->min_l :p->min_l;
		s->max_l = s->max_l > p->max_l ?s->max_l :p->max_l;
		s->reads_N += p->reads_N;
		s->bases_N += p->bases_N;
		s->bases_all += p->bases_all;
		s->single_min_N = s->single_min_N < p->single_min_N ?s->single_min_N :p->single_min_N;
		s->single_max_N = s->single_max_N > p->single_max_N ?s->single_max_N :p->single_max_N;
	}
	free(w.infos);
}

/*******************
 * quality control *
 *******************/
typedef struct {
	bseq1_t *s1, *s2;
	int flt_N;
	int *flt_n;
} qc_worker_t;

static void quality_control_worker(void *data, long seq_id, int t_id) {
	qc_worker_t *w = (qc_worker_t*)data;
	bseq1_t *s1 = &w->s1[seq_id];
	bseq1_t *s2 = (w->s2 == NULL) ?NULL :&w->s2[seq_id];
	int i, N1 = 0, N2 = 0;
	for(i = 0; i < s1->l_seq; ++i) {
		N1 += (s1->seq[i] == 'N' ?1 :0);
		if(s2) N2 += (s2->seq[i] == 'N' ?1 :0);
	}
	if(N1 >= w->flt_N || N2 >= w->flt_N) {
		++w->flt_n[t_id];
		free(s1->name); free(s1->comment); free(s1->seq); free(s1->qual);
		memset(s1, 0, sizeof(*s1));
		if(s2) {
			free(s2->name); free(s2->comment); free(s2->seq); free(s2->qual);
			memset(s2, 0, sizeof(*s2));
		}
	}
}

static int quality_control(int n, bseq1_t *s1, bseq1_t *s2, int n_threads, int flt_N) {
	int i, ret = 0;
	qc_worker_t w;
	w.s1 = s1;
	w.s2 = s2;
	w.flt_N = flt_N;
	w.flt_n = calloc(n_threads, sizeof(int));
	kt_for(n_threads, quality_control_worker, &w, n);
	for(i = 0; i < n_threads; ++i) {
		ret += w.flt_n[i];
	}
	free(w.flt_n);
	return ret;
}

static void *process(void *shared, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	opt_t *opt = aux->opt;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if(step == 0) {
		ktp_data_t *ret;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs1 = bseq_read1(opt->chunk_size * opt->n_threads, &ret->n_seqs1, &aux->has_read_bytes, aux->ks1);
		if(aux->ks2) {
			ret->seqs2 = bseq_read1(opt->chunk_size * opt->n_threads, &ret->n_seqs2, &aux->has_read_bytes, aux->ks2);
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
		aux->n_processed += data->n_seqs1;
		if(opt->qual_ctl) {
			aux->flt_n += quality_control(data->n_seqs1, data->seqs1, data->seqs2, opt->n_threads, opt->flt_N);
		}
		seqs_info(&aux->seqs1_info, opt->n_threads, data->n_seqs1, data->seqs1);
		if(aux->ks2) {
			seqs_info(&aux->seqs2_info, opt->n_threads, data->n_seqs2, data->seqs2);
		}
		aux->has_read_bytes = aux->has_read_bytes <= 0 ?(aux->f1_size + aux->f2_size) :aux->has_read_bytes;
		progress_show(&aux->bar, 1.0 * aux->has_read_bytes / (aux->f1_size + aux->f2_size));
		if(aux->has_read_bytes == aux->f1_size + aux->f2_size) {
			fprintf(stderr, "\n");
		}
//		fprintf(stderr, "[%s] processed %d sequences, %ld %% of the input file\n", __func__, data->n_seqs1 + data->n_seqs2,
//				100 * aux->has_read_bytes / (aux->f1_size + aux->f2_size));
		return data;
	} else if(step == 2) {
		if(opt->qual_ctl) {
			for(i = 0; i < data->n_seqs1; ++i) {
				const bseq1_t *p = &data->seqs1[i];
				if(p->name) {
					fprintf(opt->output_f1, "@%s/1\n", p->name);
					fprintf(opt->output_f1, "%s\n", p->seq);
					fprintf(opt->output_f1, "+\n");
					fprintf(opt->output_f1, "%s\n", p->qual);
				}
			}
			if(aux->ks2) {
				for(i = 0; i < data->n_seqs2; ++i) {
					const bseq1_t *p = &data->seqs2[i];
					if(p->name) {
						fprintf(opt->output_f2, "@%s/1\n", p->name);
						fprintf(opt->output_f2, "%s\n", p->seq);
						fprintf(opt->output_f2, "+\n");
						fprintf(opt->output_f2, "%s\n", p->qual);
					}
				}
			}
		}
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
static void output(ktp_aux_t *a) {
	if(a->opt->qual_ctl) {
		fprintf(stderr, "has filtered out %d reads of N >= %d\n", a->flt_n, a->opt->flt_N);
	}
	output_seqs_info(&a->seqs1_info);
	if(a->ks2) output_seqs_info(&a->seqs1_info);
}

int hfastq_main(int argc, char **argv) {
	if(argc == 1) {
		return usage();
	}
	int c, lo_index;
	const char *short_opts = "l:r:K:t:";
	const struct option long_opts[] = {
			// {name, has_arg, flag, val}
			{"i1", required_argument, NULL, 0},
			{"i2", required_argument, NULL, 0},
			{"o1", required_argument, NULL, 0},
			{"o2", required_argument, NULL, 0},
			{"io1", 0, NULL, 0},
			{"qc", required_argument, NULL, 0},
			{NULL, 0, NULL, 0}
	};
	opt_t *opt = opt_init();
	ktp_aux_t aux; aux_init(&aux);
	int fd1 = 0; gzFile fp1 = 0; void *ko1 = NULL;
	int fd2 = 0; gzFile fp2 = 0; void *ko2 = NULL;
	int l = -1, r = -1, no_mul_io = 0;
	while((c=getopt_long(argc, argv, short_opts, long_opts, &lo_index)) >= 0) {
		switch (c) {
			case 0:
				if(!strcmp(long_opts[lo_index].name, "i1")) {
					ko1 = kopen(optarg, &fd1);
				} else if(!strcmp(long_opts[lo_index].name, "i2")) {
					ko2 = kopen(optarg, &fd2);
				} else if(!strcmp(long_opts[lo_index].name, "o1")) {
					opt->output_f1 = fopen(optarg, "w");
				} else if(!strcmp(long_opts[lo_index].name, "o2")) {
					opt->output_f2 = fopen(optarg, "w");
				} else if(!strcmp(long_opts[lo_index].name, "io1")) {
					no_mul_io = 1;
				} else if(!strcmp(long_opts[lo_index].name, "qc")) {
					opt->qual_ctl = 1;
					opt->flt_N = atoi(optarg);
				}
				break;
			case 'l': l = atoi(optarg); break;
			case 'r': r = atoi(optarg); break;
			case 'K': opt->chunk_size = atoi(optarg); break;
			case 't': opt->n_threads = atoi(optarg); break;
			case '?': return usage();
			default : return usage();
		}
	}
	if(ko1 == NULL) {
		fprintf(stderr, "at least one FASTA/Q file\n");
		return 1;
	}
	fp1 = gzdopen(fd1, "r"); aux.ks1 = kseq_init(fp1); aux.f1_size = get_file_size(fd1);
	if(ko2 != NULL) {
		fp2 = gzdopen(fd2, "r"); aux.ks2 = kseq_init(fp2); aux.f2_size = get_file_size(fd2);
	}
	aux.opt = opt;

	kt_pipeline(no_mul_io == 1 ?1 :2, process, &aux, 3);
	output(&aux);

	aux_detroy(&aux);
	kseq_destroy(aux.ks1); gzclose(fp1); kclose(ko1);
	if(ko2 != NULL) {
		kseq_destroy(aux.ks2); gzclose(fp2); kclose(ko2);
	}
	free(opt);
	return 0;
}