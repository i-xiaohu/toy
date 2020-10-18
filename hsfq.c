//
// Created by 63175 on 2019/11/23.
//

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <zlib.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>

#include "kseq.h"
#include "sync_pe.h"
#include "ksort.h"
#include "utils.h"
#include "table.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

KSEQ_DECLARE(gzFile)
void *kopen(const char *fn, int *_fd);
int kclose(void *a);

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    hsfq         handle reodered FASTQ file\n");
	fprintf(stderr, "Usage:      hsfq         [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -i STR       input FASTQ filename\n");
	fprintf(stderr, "            -p           my smart paired-end file\n");
	fprintf(stderr, "            --check      check PE reads in fixed buffer\n");
	fprintf(stderr, "            --harc       build CS with harc's majority rule.\n");
	fprintf(stderr, "            --harc2      build Ref harc's hash principle.\n");
	fprintf(stderr, "            -K INT       pairing buffer size [10000000], only in PE mode\n");
	fprintf(stderr, "\n");
	return 1;
}

typedef struct {
	int l, r; // [l, r)
	int n; // consensus sequence length
	uint8_t *s; // consensus sequence generated form the reads batch
} cs_t;
typedef kvec_t(cs_t) cs_v;

typedef struct {
	int seqs;
	int seqs_in_cs; // in long cs
	int short_cs;
	int cs[6];
	int long_cs;
	long cs_len;
	long cs_bases;
	long all_bases;
	long dis;
} stat_t;
static stat_t hsfq_stat;

static inline void reverse(char *des, const char *src, int len) {
	int i;
	for(i = 0; i < len; ++i) {
		if(src[i] == 'A') des[len-1-i] = 'T';
		else if(src[i] == 'C') des[len-1-i] = 'G';
		else if(src[i] == 'G') des[len-1-i] = 'C';
		else des[len-1-i] = 'A';
	}
}

static void consensus_seqs(cs_v *v, int n, bseq1_t *seqs, int *dis) {
	hsfq_stat.seqs += n;
	int i, j, d, k, max_errors = 8, min_matches = 20;
	cs_t tmp; memset(&tmp, 0, sizeof(tmp));
	for(i = 0; i < n; ++i) {
		dis[i] = INT_MAX;
		for(j = 0; j < seqs[i].l_seq; ++j) {
			seqs[i].seq[j] = (seqs[i].seq[j] == 'N') ?((char)('A')) :seqs[i].seq[j];
		}
		if(i == 0) continue;
		int min_l = (seqs[i].l_seq < seqs[i-1].l_seq) ?seqs[i].l_seq :seqs[i-1].l_seq;
		for(d = 0; d < min_l - min_matches; ++d) { // enumerate distance
			int e = 0, m = 0;
			for(k = 0; k < seqs[i].l_seq && d + k < seqs[i-1].l_seq; ++k) { // enumerate bases
				if(seqs[i].seq[k] != seqs[i-1].seq[d + k]) {
					++e;
					if(e > max_errors) {
						break;
					}
				} else {
					++m;
				}
			}
			if(e <= max_errors && m >= min_matches) {
				dis[i] = d;
				break;
			}
		}
		// fixme: strand-correction should be noted.
		// if(dis[i-1] == INT_MAX), maybe it could be saved by read(i)
		if(dis[i] == INT_MAX) {
			// Job Save, to let the consensus sequence is long enough,
			// contiguous reads may share the largest overlap with too much mismatches,
			// so we try to connect with more previous reads with secondary largest overlap.
			// HARC's majority principle might work as well.
			int min_d = 0;
			for(j = i - 2; j >= 0; --j) {
				if(dis[j + 1] == INT_MAX) { break; } // read(j + 1) is not connected with read(j)
				// try to connect read(i) and read(j) to infer dis(i-1, i)
				min_d += dis[j + 1]; // min_d = dis(i-1, j)
				min_l = (seqs[i].l_seq < seqs[j].l_seq) ?seqs[i].l_seq :seqs[j].l_seq;
				if(min_d >= min_l - min_matches) { break; } // distance between read(j) and read(i) is too far.
				for(d = min_d; d < min_l - min_matches; ++d) {
					int e = 0, m = 0;
					for(k = 0; k < seqs[i].l_seq && d + k < seqs[j].l_seq; ++k) { // enumerate bases
						if(seqs[i].seq[k] != seqs[j].seq[d + k]) {
							++e;
							if(e > max_errors) {
								break;
							}
						} else {
							++m;
						}
					}
					if(e <= max_errors && m >= min_matches) {
						dis[i] = d - min_d; // d is dis(i, j), so dis(i-1, i) = dis(i, j) - dis(i-1, j)
						break;
					}
				}
				if(dis[i] != INT_MAX) {
					break;
				}
			}
		}
		if(dis[i] == INT_MAX) {
			tmp.r = i;
			if(tmp.r - tmp.l >= 5) {
				kv_push(cs_t, *v, tmp);
				++hsfq_stat.long_cs;
				hsfq_stat.seqs_in_cs += tmp.r - tmp.l;
			} else {
				++hsfq_stat.short_cs;
			}
			tmp.l = i;
		}
	}
	tmp.r = n;
	if(tmp.r - tmp.l >= 5) {
		kv_push(cs_t, *v, tmp);
		++hsfq_stat.long_cs;
		hsfq_stat.seqs_in_cs += tmp.r - tmp.l;
	} else {
		++hsfq_stat.short_cs;
	}
}

static inline uint8_t b2n(char base) {
	if(base == 'A' || base == 'a') return 0;
	else if(base == 'C' || base == 'c') return 1;
	else if(base == 'G' || base == 'g') return 2;
	else return 3;
}

static inline char n2b(uint8_t num) {
	if(num == 0) return 'A';
	else if(num == 1) return 'C';
	else if(num == 2) return 'G';
	else return 'T';
}

static void consensus_seqs_harc(cs_v *v, int n, bseq1_t *seqs, int *dis) {
	hsfq_stat.seqs += n;
	int i, j, d, k, max_errors = 8, min_matches = 20;
	cs_t tmp; memset(&tmp, 0, sizeof(tmp));
	int *cnt[4], ref_l = seqs[0].l_seq; // reads length should be equal.
	char *ref = malloc(ref_l * sizeof(char));
	char *rev = malloc(ref_l * sizeof(char));
	for(i = 0; i < 4; ++i) {
		cnt[i] = malloc(ref_l * sizeof(int));
	}
	for(i = 0; i < n; ++i) {
		dis[i] = INT_MAX;
		for(j = 0; j < seqs[i].l_seq; ++j) {
			seqs[i].seq[j] = (seqs[i].seq[j] == 'N') ?((char)('A')) :seqs[i].seq[j];
		}
		if(i > 0) {
			// try to connect to ref
			for(d = 0; d <= ref_l >> 1; ++d) { // enumerate distance
				int e = 0, m = 0;
				for(k = 0; k < seqs[i].l_seq && d + k < ref_l; ++k) { // enumerate bases
					if(seqs[i].seq[k] != ref[d + k]) {
						++e;
						if(e > max_errors) {
							break;
						}
					} else {
						++m;
					}
				}
				if(e <= max_errors && m >= min_matches) {
					dis[i] = d;
					break;
				}
			}
			if(dis[i] == INT_MAX) {
				reverse(rev, seqs[i].seq, seqs[i].l_seq);
				memcpy(seqs[i].seq, rev, seqs[i].l_seq * sizeof(char));
				// try to connect to ref
				for(d = 0; d <= ref_l >> 1; ++d) { // enumerate distance
					int e = 0, m = 0;
					for(k = 0; k < seqs[i].l_seq && d + k < ref_l; ++k) { // enumerate bases
						if(seqs[i].seq[k] != ref[d + k]) {
							++e;
							if(e > max_errors) {
								break;
							}
						} else {
							++m;
						}
					}
					if(e <= max_errors && m >= min_matches) {
						dis[i] = d;
						break;
					}
				}
				if(dis[i] == INT_MAX) {
					reverse(rev, seqs[i].seq, seqs[i].l_seq);
					memcpy(seqs[i].seq, rev, seqs[i].l_seq * sizeof(char));
				}
			}
		}
		if(dis[i] == INT_MAX) {
			// reset ref
			for(j = 0; j < 4; ++j) {
				memset(cnt[j], 0, ref_l * sizeof(int));
			}
			if(i > 0) {
				tmp.r = i;
				if(tmp.r - tmp.l >= 5) {
					kv_push(cs_t, *v, tmp);
					++hsfq_stat.long_cs;
					hsfq_stat.cs[5] += tmp.r - tmp.l;
					hsfq_stat.seqs_in_cs += tmp.r - tmp.l;
				} else {
					++hsfq_stat.short_cs;
					hsfq_stat.cs[tmp.r-tmp.l] += tmp.r - tmp.l;
				}
				tmp.l = i;
			}
		} else {
			// delete dis[i] bases front
			for(j = 0; j < ref_l - dis[i]; ++j) {
				ref[j] = ref[j + dis[i]];
				for(k = 0; k < 4; ++k) {
					cnt[k][j] = cnt[k][j + dis[i]];
				}
			}
			for(j = ref_l - dis[i]; j < ref_l; ++j) {
				for(k = 0; k < 4; ++k) {
					cnt[k][j] = 0;
				}
			}
		}
		// vote for read(i)
		for(j = 0; j < ref_l; ++j) {
			++cnt[b2n(seqs[i].seq[j])][j];
		}
		for(j = 0; j < ref_l; ++j) {
			int max = -1, max_id = 0;
			for(k = 0; k < 4; ++k) {
				if(cnt[k][j] > max) {
					max = cnt[k][j];
					max_id = k;
				}
			}
			ref[j] = n2b(max_id);
		}
	}
	tmp.r = n;
	if(tmp.r - tmp.l >= 5) {
		kv_push(cs_t, *v, tmp);
		++hsfq_stat.long_cs;
		hsfq_stat.cs[5] += tmp.r - tmp.l;
		hsfq_stat.seqs_in_cs += tmp.r - tmp.l;
	} else {
		++hsfq_stat.short_cs;
		hsfq_stat.cs[tmp.r-tmp.l] += tmp.r - tmp.l;
	}
	free(ref); free(rev);
	for(i = 0; i < 4; ++i) {
		free(cnt[i]);
	}
}

#define MER_t uint8_t
#define MER_K 4

static void ref_harc2(cs_v *v, int n, bseq1_t *seqs, int *dis) {
	/* @ follow the principle of HARC */
	/* @ memory usage of Ref <= sum of read length */
	/* @ slide K-mer forward on Ref to match index on next read */
	/* @ backward search for Rev, and Rev is after next read */
	hsfq_stat.seqs += n;
	int i, j, k, max_l = 0;
	cs_t tmp; memset(&tmp, 0, sizeof(tmp));
	for(i = 0; i < n; ++i) {
		char *seq = seqs[i].seq;
		for(j = 0; j < seqs[i].l_seq; ++j) {
			seq[j] = seq[j] == 'N' ?'A' :seq[j];
			seq[j] = seq[j] < 4 ?seq[j] :b2n(seq[j]);
		}
		max_l += seqs[i].l_seq;
	}
	uint8_t *ref = malloc(max_l * sizeof(uint8_t));
	uint16_t *cnt[4];
	for(i = 0; i < 4; ++i) {
		cnt[i] = calloc(max_l, sizeof(uint16_t));
	}
	bseq1_t *first = &seqs[0];
	for(i = 0; i < first->l_seq; ++i) {
		++cnt[first->seq[i]][i];
		ref[i] = first->seq[i];
	}
	int fix_len = first->l_seq;

	// fixme: only support reads of length >= 64
	// todo: support variable length
	tmp.l = 0; dis[0] = INT_MAX;
	int fan = 0, zheng = 0;
	uint8_t *rev = malloc(fix_len * sizeof(uint8_t));
	for(i = 1; i < n; ++i) {
		const char *seq = seqs[i].seq;
		const int l_seq = seqs[i].l_seq; assert(l_seq >= 64); assert(l_seq == 150);
		int m = l_seq >> 1;
		MER_t H1 = 0, H2 = 0; // H1=[m-32, m), H2=[m, m+32)
		for(j = 0; j < MER_K; ++j) {
			H1 <<= 2; H1 += seq[m-MER_K+j];
			H2 <<= 2; H2 += seq[m+j];
		}

		m = fix_len >> 1;
		for(j = 0; j < fix_len; ++j) {
			rev[j] = 3 - ref[fix_len-1-j];
		}
		MER_t h1 = 0, h2 = 0, r1 = 0, r2 = 0;
		for(j = 0; j < MER_K; ++j) {
			h1 <<= 2; h1 += ref[m-MER_K+j];
			h2 <<= 2; h2 += ref[m+j];
			r1 <<= 2; r1 += rev[m-MER_K+j];
			r2 <<= 2; r2 += rev[m+j];
		}
		int d = -1, rc = 0;
		for(j = 0; ; ++j) {
			if(h1 == H1 || h2 == H2) {
				d = j;
				int mis = 0;
				for(k = 0; k + d < fix_len; ++k) {
					mis += (ref[k + d] != seq[k]) ?1 :0;
				}
				if(mis <= 8) {
					++zheng;
					break;
				}
				d = -1;
			}
			if(r1 == H1 || r2 == H2) {
				d = j;
				int mis = 0;
				for(k = 0; k + d < fix_len; ++k) {
					mis += (rev[k] != seq[k + d]) ?1 :0;
				}
				if(mis <= 8) {
					++fan;
					rc = 1;
					break;
				}
				d = -1;
			}
			// forward shift
			int brk = 1;
			if(m + j < fix_len) { h1 <<= 2; h1 += ref[m + j]; brk = 0; }
			if(m + MER_K + j < fix_len) { h2 <<= 2; h2 += ref[m + MER_K + j]; brk = 0; }
			// backward shift
			if(m - MER_K - 1 - j >= 0) { r1 >>= 2; r1 += 1UL * rev[m - MER_K - 1 - j] << (MER_K*2 - 2); brk = 0; }
			if(m - 1 - j >= 0) { r2 >>= 2; r2 += 1UL * rev[m - 1 - j] << (MER_K*2 - 2); brk = 0; }
			if(brk) break;
		}
		if(d != -1) {
			dis[i] = d;
			// delete dis[i] bases front
			for(j = 0; j < fix_len - dis[i]; ++j) {
				ref[j] = ref[j + dis[i]];
				for(k = 0; k < 4; ++k) {
					cnt[k][j] = cnt[k][j + dis[i]];
				}
			}
			for(j = fix_len - dis[i]; j < fix_len; ++j) {
				for(k = 0; k < 4; ++k) {
					cnt[k][j] = 0;
				}
			}
		} else {
			// reset ref
			for(j = 0; j < 4; ++j) {
				memset(cnt[j], 0, fix_len * sizeof(int16_t));
			}
			dis[i] = INT_MAX;
			tmp.r = i;
			if(tmp.r - tmp.l >= 5) {
				kv_push(cs_t, *v, tmp);
				hsfq_stat.cs[5] += tmp.r - tmp.l;
			} else {
				hsfq_stat.cs[tmp.r-tmp.l] += tmp.r - tmp.l;
			}
			tmp.l = i;
		}
		// vote for read(i)
		for(j = 0; j < fix_len; ++j) {
			if(!rc) ++cnt[seq[j]][j];
			else ++cnt[3 - seq[fix_len-1-j]][j];
		}
		for(j = 0; j < fix_len; ++j) {
			int max = -1, max_id = 0;
			for(k = 0; k < 4; ++k) {
				if(cnt[k][j] > max) {
					max = cnt[k][j];
					max_id = k;
				}
			}
			ref[j] = max_id;
		}
	}
	tmp.r = n;
	if(tmp.r - tmp.l >= 5) {
		kv_push(cs_t, *v, tmp);
		hsfq_stat.cs[5] += tmp.r - tmp.l;
	} else {
		hsfq_stat.cs[tmp.r-tmp.l] += tmp.r - tmp.l;
	}
	free(ref); free(rev);
	for(i = 0; i < 4; ++i) {
		free(cnt[i]);
	}
}

typedef struct {
	MER_t h;
	int id;
} h_id_t;

#define h_id_lt(a, b) ((a).h < (b).h)
KSORT_INIT(HARC, h_id_t, h_id_lt)

static int inline lower_bound(int n, h_id_t *a, MER_t x) {
	int ret = n, l = 0, r = n - 1, mid;
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

#define bseq1_lt(a, b) (strcmp((a).name, (b).name) < 0)
KSORT_INIT(bseq1, bseq1_t, bseq1_lt)

static int check_name(int n, bseq1_t *seqs) {
	int i;
	for(i = 0; i < n; ++i) {
		if((i&1) == 0 && i > 0) {
			if(strcmp(seqs[i].name, seqs[i-1].name) == 0) {
				return 0;
			}
		}
		if((i&1) != 0) {
			if(strcmp(seqs[i].name, seqs[i-1].name) != 0) {
				return 0;
			}
		}
	}
	return 1;
}

static void get_coverage(cs_v *v, bseq1_t *seqs, int *dis) {
	int i ,j;
	for(i = 0; i < v->n; ++i) {
		cs_t *cs = &v->a[i];
		int cs_len = 0;
		for(j = cs->l; j < cs->r; ++j) {
			if(j == cs->l) dis[j] = 0;
			else {
				hsfq_stat.dis += dis[j];
				dis[j] = dis[j - 1] + dis[j];
			}
			bseq1_t *s = &seqs[j];
			cs_len = cs_len >= s->l_seq + dis[j] ?cs_len :s->l_seq + dis[j];
			hsfq_stat.all_bases += s->l_seq;
		}
		hsfq_stat.cs_bases += cs_len;
	}
}

static void process_pe(void *ks_, int K, int check, int c_way) {
	kseq_t *ks = (kseq_t*)ks_;
	int i, n; long bytes = 0;
	bseq1_t *seqs1 = NULL;
	memset(&hsfq_stat, 0, sizeof(hsfq_stat));
	cs_v v; kv_init(v);
	while((seqs1 = bseq_read1(K, &n, &bytes, ks)) != NULL) {
		int *dis = malloc(n * sizeof(int));
		v.n = 0;
		if(c_way == 0) {
//			consensus_seqs(&v, n, seqs1, dis);
		} else if(c_way == 1) {
//			consensus_seqs_harc(&v, n, seqs1, dis);
		} else if(c_way == 2) {
//			ref_harc2_pe(&v, n, seqs1, seq2, dis);
		}
		get_coverage(&v, seqs1, dis);
		if(check) {
			ks_introsort(bseq1, n, seqs1); // sort by qname
			assert(check_name(n, seqs1) == 1);
		}
		free(dis);
		for(i = 0; i < n; ++i) {
			free(seqs1[i].name); free(seqs1[i].comment); free(seqs1[i].seq); free(seqs1[i].qual);
		}
		free(seqs1);
	}
	free(v.a);

	assert((hsfq_stat.seqs & 1) == 0);
	/**
	 * coverage = all_bases / cs_bases
	 * distance = cs_length / seqs_in_cs
	 */
	tab_t *t = tab_init(3, 8, "consensus sequence statistics", NULL);
	char *a[8] = {"# seq-pairs", "# seqs-in-CS", "seq-dis", "# CS", "# short CS", "# long CS", "CS-length", "CS-coverage"};
	for(i = 1; i <= 8; ++i) {
		tab_fill_cell1(t, 1, i, "%s", a[i-1]);
	}
	tab_fill_cell1(t, 2, 1, "%d", hsfq_stat.seqs >> 1);
	tab_fill_cell1(t, 2, 2, "%d", hsfq_stat.seqs_in_cs); tab_fill_cell1(t, 3, 2, "%.3f %%", 100.0 * hsfq_stat.seqs_in_cs / hsfq_stat.seqs);
	tab_fill_cell1(t, 2, 3, "%.3f", 1.0 * hsfq_stat.dis / (hsfq_stat.seqs_in_cs - hsfq_stat.long_cs));
	tab_fill_cell1(t, 2, 4, "%d", hsfq_stat.short_cs + hsfq_stat.long_cs);
	tab_fill_cell1(t, 2, 5, "%d", hsfq_stat.short_cs); tab_fill_cell1(t, 3, 5, "%.3f %%", 100.0 * hsfq_stat.short_cs / (hsfq_stat.short_cs + hsfq_stat.long_cs));
	tab_fill_cell1(t, 2, 6, "%d", hsfq_stat.long_cs); tab_fill_cell1(t, 3, 6, "%.3f %%", 100.0 * hsfq_stat.long_cs / (hsfq_stat.short_cs + hsfq_stat.long_cs));
	tab_fill_cell1(t, 2, 7, "%ld", hsfq_stat.cs_bases); tab_fill_cell1(t, 3, 7, "%.3f", 1.0 * hsfq_stat.cs_bases / hsfq_stat.long_cs);
	tab_fill_cell1(t, 2, 8, "%.3f X", 1.0 * hsfq_stat.all_bases / hsfq_stat.cs_bases);
	tab_display(t);
	tab_destroy(t);
	free(t);
}

static void process_se(void *ks_, int K, int c_way) {
	kseq_t *ks = (kseq_t*)ks_;
	int i, n; long bytes = 0;
	bseq1_t *seqs = NULL;
	memset(&hsfq_stat, 0, sizeof(hsfq_stat));
	cs_v v; kv_init(v);
	while((seqs = bseq_read1(K, &n, &bytes, ks)) != NULL) {
		int *dis = malloc(n * sizeof(int));
		v.n = 0;
		if(c_way == 0) {
			consensus_seqs(&v, n, seqs, dis);
		} else if(c_way == 1) {
			consensus_seqs_harc(&v, n, seqs, dis);
		} else if(c_way == 2) {
			ref_harc2(&v, n, seqs, dis);
		}
		get_coverage(&v, seqs, dis);
		free(dis);
		for(i = 0; i < n; ++i) {
			free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual);
		}
		free(seqs);
	}
	free(v.a);
	int r = 3, c = 5;
	char *hdr[25] = {"1", "2", "3", "4", ">=5"};
	tab_t *t = tab_init(r, c, "CS length", NULL);
	int cs_n = 0;
	for(i = 1; i <= c; ++i) {
		tab_fill_cell1(t, 1, i, "%s", hdr[i-1]);
		tab_fill_cell1(t, 2, i, "%d", hsfq_stat.cs[i]);
		tab_fill_cell1(t, 3, i, "%.2f %%", 100.0 * hsfq_stat.cs[i] / hsfq_stat.seqs);
		cs_n += hsfq_stat.cs[i];
	}
	assert(cs_n == hsfq_stat.seqs);
	tab_display(t); tab_destroy(t);

	t = tab_init(3, 8, "consensus sequence statistics", NULL);
	char *a[8] = {"# seqs", "# seqs-in-CS", "seq-dis", "# CS", "# short CS", "# long CS", "CS-length", "CS-coverage"};
	for(i = 1; i <= 8; ++i) {
		tab_fill_cell1(t, 1, i, "%s", a[i-1]);
	}
	tab_fill_cell1(t, 2, 1, "%d", hsfq_stat.seqs);
	tab_fill_cell1(t, 2, 2, "%d", hsfq_stat.seqs_in_cs); tab_fill_cell1(t, 3, 2, "%.3f %%", 100.0 * hsfq_stat.seqs_in_cs / hsfq_stat.seqs);
	tab_fill_cell1(t, 2, 3, "%.3f", 1.0 * hsfq_stat.dis / (hsfq_stat.seqs_in_cs - hsfq_stat.long_cs));
	tab_fill_cell1(t, 2, 4, "%d", hsfq_stat.short_cs + hsfq_stat.long_cs);
	tab_fill_cell1(t, 2, 5, "%d", hsfq_stat.short_cs); tab_fill_cell1(t, 3, 5, "%.3f %%", 100.0 * hsfq_stat.short_cs / (hsfq_stat.short_cs + hsfq_stat.long_cs));
	tab_fill_cell1(t, 2, 6, "%d", hsfq_stat.long_cs); tab_fill_cell1(t, 3, 6, "%.3f %%", 100.0 * hsfq_stat.long_cs / (hsfq_stat.short_cs + hsfq_stat.long_cs));
	tab_fill_cell1(t, 2, 7, "%ld", hsfq_stat.cs_bases); tab_fill_cell1(t, 3, 7, "%.3f", 1.0 * hsfq_stat.cs_bases / hsfq_stat.long_cs);
	tab_fill_cell1(t, 2, 8, "%.3f X", 1.0 * hsfq_stat.all_bases / hsfq_stat.cs_bases);
	tab_display(t); tab_destroy(t);
	free(t);
}

int hsfq_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c, lo_index;
	const char *short_opts = "i:K:p";
	const struct option long_opts[] = {
			// {name, has_arg, flag, val}
			{"check", 0, NULL, 0},
			{"harc", 0, NULL, 0},
			{"harc2", 0, NULL, 0},
			{NULL, 0, NULL, 0}
	};
	int pe = 0, K = 10000000, check = 0, connect_way = 0;
	int fd = 0; gzFile fp = 0; void *ko = NULL;
	while((c=getopt_long(argc, argv, short_opts, long_opts, &lo_index)) >= 0) {
		switch (c) {
			case 0:
				if(!strcmp(long_opts[lo_index].name, "check")) {
					check = 1;
				} else if(!strcmp(long_opts[lo_index].name, "harc")) {
					connect_way = 1;
				} else if(!strcmp(long_opts[lo_index].name, "harc2")) {
					connect_way = 2;
				}
				break;
			case 'i': ko = kopen(optarg, &fd); break;
			case 'K': K = atoi(optarg); break;
			case 'p': pe = 1; break;
			case '?': return usage();
			default : return usage();
		}
	}
	if(ko == NULL) {
		fprintf(stderr, "can't open file\n");
		return 1;
	}
	fp = gzdopen(fd, "r"); kseq_t *ks = kseq_init(fp);
	if(pe) {
		process_pe(ks, K, check, connect_way);
	} else {
		process_se(ks, K, connect_way);
	}
	kseq_destroy(ks); gzclose(fp); kclose(ko);
	return 0;
}
