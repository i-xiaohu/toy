//
// Created by 63175 on 2019/12/14.
//

#include <stdio.h>
#include <zlib.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>

#include "kseq.h"
#include "utils.h"
#include "progress.h"
#include "sync_pe.h"
#include "kvec.h"
#include "table.h"
#include "ksort.h"

#define MATCH_COEF 0.10
#define ERROR_COEF 0.08
#define KMER_L     32

KSEQ_DECLARE(gzFile)
extern void *kopen(const char *fn, int *_fd);
extern int kclose(void *a);
extern unsigned char nst_nt4_table[256];

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

/*
 * todo: use real hash instead of binary search
 * todo: implement pair-end read connect
 */

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    reorder      test reorder performance of Spring(reorder).\n");
	fprintf(stderr, "Usage:      reorder      [options] <in1.fq> [in2.fq]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            --force      enumerate distance between consecutive reads.\n");
	fprintf(stderr, "            --hash       build CS with K-Mer hash.\n");
	fprintf(stderr, "            --mtry       max search reads number\n");
	fprintf(stderr, "            --io1        no multiple IO thread in pipeline\n");
	fprintf(stderr, "            -K INT       buffer size [10000000]\n");
	fprintf(stderr, "            -t INT       threads [8]\n");
	fprintf(stderr, "\n");
	return 1;
}

typedef struct {
	int way;
	int chunk_size;
	int mtry;
	int n_threads;
} opt_t;
static opt_t *opt_init() {
	opt_t *o = calloc(1, sizeof(opt_t));
	o->way = 0;
	o->chunk_size = 10000000;
	o->mtry = 8;
	o->n_threads = 8;
	return o;
}

typedef struct {
	opt_t *opt;
	kseq_t *ks1, *ks2;
	long f1_size, f2_size, has_read_bytes;
	long n_processed;
	progress_t *bar;
	int cnt[6];
} ktp_aux_t;

typedef struct {
	int n_seqs1, n_seqs2;
	bseq1_t *seqs1, *seqs2;
} ktp_data_t;

typedef struct {
	int l, r;
	size_t bases;
} cs_t;
typedef kvec_t(cs_t) cs_v;

typedef struct {
	int *stop;  // thread i processes [stop[i-1], stop[i])
	bseq1_t *seqs1, *seqs2;
	cs_v *sub_cs, cs;
	int **sub_cnt; // consensus length statistics
	int *offset; // offset in consensus sequence
	int max_try;
} worker_t;

/** build ref to connect next read */
typedef struct {
	int n, m;
	uint8_t *a;
	int *cnt[4];
} ref_t;

static ref_t *ref_init() {
	ref_t *r = calloc(1, sizeof(ref_t));
	return r;
}

static void ref_free(ref_t *r) {
	free(r->a);
	free(r->cnt[0]); free(r->cnt[1]); free(r->cnt[2]); free(r->cnt[3]);
	free(r);
}

static inline void ref_add(ref_t *r, int d, int l, const uint8_t *seq, int rev) {
	assert(d <= r->n);
	int i, n;
	n = (d + l > r->n) ?(d + l) :r->n;
	if(n >= r->m) {
		r->m = n << 1;
		r->a = realloc(r->a, r->m * sizeof(uint8_t));
		memset(r->a + r->n, 0, (r->m - r->n) * sizeof(uint8_t));
		for(i = 0; i < 4; ++i) {
			r->cnt[i] = realloc(r->cnt[i], r->m * sizeof(int));
			memset(r->cnt[i] + r->n, 0, (r->m - r->n) * sizeof(int));
		}
	}
	r->n = n;
	for(i = 0; i < l; ++i) {
		uint8_t c = rev == 1 ?3-seq[l-1-i] :seq[i];
		++r->cnt[c][d+i];
		if(r->cnt[c][d+i] > r->cnt[r->a[d+i]][d+i]) {
			r->a[d+i] = c;
		}
	}
}

static inline void ref_delete(ref_t *r, int d) {
	assert(d <= r->n);
	int i;
	for(i = 0; i < r->n - d; ++i) {
		r->a[i] = r->a[i+d];
		r->cnt[0][i] = r->cnt[0][i+d];
		r->cnt[1][i] = r->cnt[1][i+d];
		r->cnt[2][i] = r->cnt[2][i+d];
		r->cnt[3][i] = r->cnt[3][i+d];
	}
	r->n -= d;
	memset(r->a + r->n, 0, (r->m - r->n) * sizeof(uint8_t));
	for(i = 0; i < 4; ++i) {
		memset(r->cnt[i] + r->n, 0, (r->m - r->n) * sizeof(int));
	}
}

/** build KMer Hash table for query seq */
typedef struct {
	uint64_t h;
	int p;
} kh_t;
typedef kvec_t(kh_t) kh_v;

#define kh_lt(a, b) ((a).h != (b).h ?(a).h < (b).h :(a).p < (b).p)
KSORT_INIT(kh, kh_t, kh_lt)

static inline void kh_seq(int l, const uint8_t *seq, kh_v *v) {
	v->n = 0;
	int i, j;
	for(i = 0; i < l - KMER_L; ++i) {
		kh_t t; t.p = i; t.h = 0;
		for(j = 0; j < KMER_L; ++j) {
			t.h <<= 2;
			t.h += seq[i + j];
		}
		kv_push(kh_t, *v, t);
	}
	ks_introsort(kh, v->n, v->a);
}

static inline int binary_search(int n, const kh_t *a, uint64_t x) {
	int l = 0, r = n, mid, ret = n;
	while(l <= r) {
		mid = (l + r) >> 1;
		if(a[mid].h >= x) {
			ret = mid;
			r = mid - 1;
		} else {
			l = mid + 1;
		}
	}
	return ret;
}

static inline int kh_match(int l1, const uint8_t *s1, int l2, const uint8_t *s2, int d) {
	int min_l = (l1 < l2) ?l1 :l2;
	int max_errors = (int)(ERROR_COEF * min_l + 0.499);
	int min_matches = (int)(MATCH_COEF * min_l + 0.499); min_matches = min_matches < 19 ?19 :min_matches;
	int k, e = 0, m = 0;
	for(k = 0; k < l2 && k + d < l1; ++k) {
		if(s2[k] != s1[d + k]) { ++e; }
		else { ++m; }
		if(e > max_errors) { break; }
	}
	if(e <= max_errors && m >= min_matches) { return 1; }
	return 0;
}

static inline int kh_find(ref_t *ref, int len, const uint8_t *seq, kh_v *v, int max_try) {
	// time complexity:
	//   L denotes maximal length of reads
	//   D denotes average distance between consecutive reads
	//   2 * (L * log(2, L) + D * L + 2 * L)
	//   = (D + 2 + log(2, L) * L * 2
	// for reads of length 256, distance 10, time cost is 40 * L
	// in the force way, time cost is D * L, which is much faster than kmer-hash
	int i, j, k, ret = INT_MAX;
	for(k = 0; k < 2; ++k) {
		// k = 0, forward search; k = 1, reverse complement search
		for(i = 0; i < ref->n - KMER_L; ++i) {
			uint64_t x = 0;
			for(j = 0; j < KMER_L; ++j) {
				x <<= 2; x += ref->a[i + j];
			}
			int res = binary_search(v->n, v->a, x);
			for(j = res; j < v->n; ++j) {
				if(j - res + 1 > max_try) { break; }
				const kh_t *p = &v->a[j];
				if(p->h != x) { break; }
				int d = (k == 0) ?(i - p->p) :(p->p - i);
				if(d < 0) { continue; } // next read have to be after ref
				if(k == 0 && kh_match(ref->n, ref->a, len, seq, d)) { ret = d; break; }
				if(k == 1 && kh_match(len, seq, ref->n, ref->a, d)) { ret = -d; break; }
			}
			if(ret != INT_MAX) break;
		}
		if(k == 0 && ret != INT_MAX) break;
		for(j = 0; j < ref->n >> 1; ++j) {
			uint8_t swap = ref->a[j];
			ref->a[j] = ref->a[ref->n-1-j];
			ref->a[ref->n-1-j] = swap;
		}
		for(j = 0; j < ref->n; ++j) {
			ref->a[j] = (uint8_t)(3 - ref->a[j]);
		}
	}
	return ret;
}

/** connect consecutive reads by Kmer-Hash */
static void kmer_hash(void *data, long seq_id, int t_id) {
	assert(seq_id == t_id);
	worker_t *w = (worker_t*)data;
	int start = t_id == 0 ?0 :w->stop[t_id-1], end = w->stop[t_id];
	bseq1_t *seqs1 = w->seqs1;
	int *offset = w->offset;
	cs_v *v = &w->sub_cs[t_id];
	int *cnt = w->sub_cnt[t_id];

	int i, j, d;
	cs_t cs; cs.l = start; cs.r = end;
	ref_t *r = ref_init();
	kh_v hv; kv_init(hv);
	for(i = start; i < end; ++i) {
		uint8_t *seq = (uint8_t*)seqs1[i].seq;
		int l_seq = seqs1[i].l_seq;
		offset[i] = -1;
		for(j = 0; j < l_seq; ++j) {
			seq[j] = (seq[j] == 'N') ?(uint8_t)('A') :seq[j];
			seq[j] = seq[j] < 4 ?seq[j] :nst_nt4_table[seq[j]];
		}
		if(i == start) {
			ref_add(r, 0, l_seq, seq, 0);
			continue;
		}
		int is_rev = 0;
		// find all kmers on seq
		kh_seq(l_seq, seq, &hv);
		d = kh_find(r, l_seq, seq, &hv, w->max_try);
		offset[i] = d == INT_MAX ?-1 :(d < 0 ?-d :d);
		is_rev = d < 0 ?1 :0;
		if(offset[i] == -1 || i == end - 1) {
			cs.r = i == end-1 ?end :i;
			int cs_n = cs.r - cs.l;
			if(cs_n >= 5) { kv_push(cs_t, *v, cs); cnt[5] += cs_n; }
			else { cnt[cs_n] += cs_n; }
			cs.l = i;
			ref_delete(r, r->n);
			ref_add(r, 0, l_seq, seq, is_rev);
		} else {
			ref_delete(r, offset[i]);
			ref_add(r, 0, l_seq, seq, is_rev);
		}
	}
	free(hv.a);
	ref_free(r);
}

/** connect consecutive reads forcely */
static void force(void *data, long seq_id, int t_id) {
	assert(seq_id == t_id);
	worker_t *w = (worker_t*)data;
	int start = t_id == 0 ?0 :w->stop[t_id-1], end = w->stop[t_id];
	bseq1_t *seqs = w->seqs1;
	int *offset = w->offset;
	cs_v *v = &w->sub_cs[t_id];
	int *cnt = w->sub_cnt[t_id];

	int i, j, d, k;
	cs_t cs; cs.l = start; cs.r = end;
	ref_t *r = ref_init();
	for(i = start; i < end; ++i) {
		uint8_t *seq = (uint8_t*)seqs[i].seq;
		int l_seq = seqs[i].l_seq;
		offset[i] = -1;
		for(j = 0; j < l_seq; ++j) {
			seq[j] = (seq[j] == 'N') ?(uint8_t)('A') :seq[j];
			seq[j] = seq[j] < 4 ?seq[j] :nst_nt4_table[seq[j]];
		}
		if(i == start) {
			ref_add(r, 0, l_seq, seq, 0);
			continue;
		}
		int min_l = (l_seq < r->n) ?l_seq :r->n;
		int max_errors = (int)(ERROR_COEF * min_l + 0.499);
		int min_matches = (int)(MATCH_COEF * min_l + 0.499); min_matches = min_matches < 19 ?19 :min_matches;
		int is_rev = 0;
		for(d = 0; d < min_l - min_matches; ++d) { // enumerate distance
			int e1 = 0, m1 = 0, e2 = 0, m2 = 0;
			for(k = 0; k < l_seq && d + k < r->n; ++k) {
				if(seq[k] != r->a[d + k]) { ++e1; }
				else { ++m1; }
				// reverse seq
				if(3 - seq[l_seq-1-k] != r->a[d + k]) { ++e2; }
				else { ++m2; }
				if(e1 > max_errors && e2 > max_errors) { break; }
			}
			if(e1 <= max_errors && m1 >= min_matches) {
				offset[i] = d; break;
			} else if(e2 <= max_errors && m2 >= min_matches) {
				is_rev = 1; offset[i] = d; break;
			}
		}
		if(offset[i] == -1 || i == end - 1) {
			cs.r = i == end-1 ?end :i;
			int cs_n = cs.r - cs.l;
			if(cs_n >= 5) { kv_push(cs_t, *v, cs); cnt[5] += cs_n; }
			else { cnt[cs_n] += cs_n; }
			cs.l = i;
			ref_delete(r, r->n);
			ref_add(r, 0, l_seq, seq, is_rev);
		} else {
			ref_delete(r, offset[i]);
			ref_add(r, 0, l_seq, seq, is_rev);
		}
	}
	ref_free(r);
}

static void reorder(ktp_aux_t *aux, ktp_data_t *data) {
	int n = data->n_seqs1;
	bseq1_t *seqs1 = data->seqs1;
	bseq1_t *seqs2 = data->seqs2;
	opt_t *opt = aux->opt;
	int chunk_size = opt->chunk_size;
	int way = opt->way;

	int i, j, n_threads = 0;
	worker_t w; memset(&w, 0, sizeof(w));
	w.max_try = opt->mtry;
	w.stop = malloc(opt->n_threads * sizeof(int));
	int bytes = 0, cnt_reads = 0;
	for(i = 0; i < n; ++i) {
		++cnt_reads;
		bytes += seqs1[i].l_seq;
		if(i == n - 1 || (bytes >= chunk_size && (cnt_reads & 1) == 0)) {
			w.stop[n_threads++] = i + 1;
			cnt_reads = 0;
			bytes = 0;
		}
	}
	assert(n_threads <= opt->n_threads);
	w.seqs1 = seqs1; w.seqs2 = seqs2;
	w.sub_cs = calloc(n_threads, sizeof(cs_v));
	for(i = 0; i < n_threads; ++i) {
		w.sub_cs[i].m = w.stop[i] - (i == 0 ?0 :w.stop[i-1]);
		w.sub_cs[i].a = malloc(w.sub_cs[i].m * sizeof(cs_t));
	}
	w.cs.m = n; w.cs.a = malloc(n * sizeof(cs_t));
	w.sub_cnt = calloc(n_threads, sizeof(int*));
	for(i = 0; i < n_threads; ++i) { w.sub_cnt[i] = calloc(6, sizeof(int)); }
	w.offset = malloc(n * sizeof(int));

	if(way == 0) {
		kt_for(n_threads, force, &w, n_threads);
	} else if(way == 1) {
		kt_for(n_threads, kmer_hash, &w, n_threads);
	}

	// sum up subtask
	for(i = 0; i < n_threads; ++i) {
		const cs_v *p = &w.sub_cs[i];
		for(j = 0; j < p->n; ++j) {
			kv_push(cs_t, w.cs, p->a[j]);
		}
		for(j = 1; j <= 5; ++j) {
			aux->cnt[j] += w.sub_cnt[i][j];
		}
	}

	// free allocated memory
	free(w.offset);
	for(i = 0; i < n_threads; ++i) { free(w.sub_cnt[i]); }
	free(w.sub_cnt);
	free(w.cs.a);
	for(i = 0; i < n_threads; ++i) { free(w.sub_cs[i].a); }
	free(w.sub_cs);
	free(w.stop);
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
		reorder(aux, data);
		aux->n_processed += data->n_seqs1;
		aux->has_read_bytes = aux->has_read_bytes <= 0 ?(aux->f1_size + aux->f2_size) :aux->has_read_bytes;
		progress_show(aux->bar, 1.0 * aux->has_read_bytes / (aux->f1_size + aux->f2_size));
		if(aux->has_read_bytes == aux->f1_size + aux->f2_size) {
			fprintf(stderr, "\n");
		}
		return data;
	} else if(step == 2) {
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

static void output(ktp_aux_t *a) {
	fprintf(stderr, "[%s] processed %ld sequences\n", __func__, a->n_processed);
	int i, n = 3, m = 5;
	tab_t *t = tab_init(n, m, "cs statistics", NULL);
	char *hdr[25] = {"= 1", "= 2", "= 3", "= 4", ">= 5"};
	for(i = 1; i <= m; ++i) {
		tab_fill_cell1(t, 1, i, "%s", hdr[i-1]);
		tab_fill_cell1(t, 2, i, "%d", a->cnt[i]);
		tab_fill_cell1(t, 3, i, "%.2f %%", 100.0 * a->cnt[i] / a->n_processed);
	}
	tab_display(t); tab_destroy(t);
}

int spring_reorder_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c, lo_index, no_mul_io = 0;
	const char *short_opts = "K:t:";
	const struct option long_opts[] = {
			{"force", 0, NULL, 0},
			{"hash", 0, NULL, 0},
			{"mtry", required_argument, NULL, 0},
			{"io1", 0, NULL, 0},
			{NULL, 0, NULL, 0}
	};
	ktp_aux_t aux; memset(&aux, 0, sizeof(aux));
	opt_t *opt;
	progress_t *bar;
	aux.opt = opt = opt_init();
	aux.bar = bar = calloc(1, sizeof(progress_t)); progress_init(bar, "", 100, PROGRESS_CHR_STYLE);
	int fd1 = 0; gzFile fp1 = 0; void *ko1 = NULL;
	int fd2 = 0; gzFile fp2 = 0; void *ko2 = NULL;
	while((c=getopt_long(argc, argv, short_opts, long_opts, &lo_index)) >= 0) {
		switch (c) {
			case 0:
				if(!strcmp(long_opts[lo_index].name, "force")) {
					opt->way = 0;
				} else if(!strcmp(long_opts[lo_index].name, "hash")) {
					opt->way = 1;
				} else if(!strcmp(long_opts[lo_index].name, "mtry")) {
					opt->mtry = atoi(optarg);
				} else if(!strcmp(long_opts[lo_index].name, "io1")) {
					no_mul_io = 1;
				}
				break;
			case 'K': opt->chunk_size = atoi(optarg); break;
			case 't': opt->n_threads = atoi(optarg); break;
			default : return usage();
		}
	}
	if(optind > argc-1 || optind + 2 <= argc-1) {
		return usage();
	}

	ko1 = kopen(argv[optind], &fd1);
	if(ko1 == NULL) { fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind]); }
	fp1 = gzdopen(fd1, "r"); aux.ks1 = kseq_init(fp1); aux.f1_size = get_file_size(fd1);
	if(optind + 1 <= argc-1) {
		ko2 = kopen(argv[optind + 1], &fd2);
		if(ko2 == NULL) { fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]); }
		fp2 = gzdopen(fd2, "r"); aux.ks2 = kseq_init(fp2); aux.f2_size = get_file_size(fd2);
	}

	kt_pipeline(no_mul_io == 1 ?1 :2, process, &aux, 3);

	output(&aux);

	// free allocated memory
	if(ko2 != NULL) {
		kseq_destroy(aux.ks2); gzclose(fp2); kclose(ko2);
	}
	kseq_destroy(aux.ks1); gzclose(fp1); kclose(ko1);
	progress_destroy(bar); free(bar);
	free(opt);
	return 0;
}
