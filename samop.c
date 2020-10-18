//
// Created by 63175 on 2019/9/19.
//

#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "samop.h"
#include "utils.h"
#include "ksort.h"
#include "kstring.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    samop handle SAM file\n");
	fprintf(stderr, "Usage:      samop [options] <in.sam>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -i STR        input SAM filename\n");
//	fprintf(stderr, "            -o STR        output SAM filename\n");
	fprintf(stderr, "            -r INT        show a (pair of) SAM record\n");
	fprintf(stderr, "            --stat        the statistics of SAM reads\n");
	fprintf(stderr, "            --wgsim-eval  \n");
	fprintf(stderr, "            --cov         the coverage of mapped reads\n");
	fprintf(stderr, "            --dis         the average distance of mapped reads\n");
//	fprintf(stderr, "            --pri         discard supplementary and secondary alignments, only keep unmapped and primary hits.\n");
//	fprintf(stderr, "            --cs          consensus sequences built from coordinate-sorted reads\n");
	fprintf(stderr, "\n");
	return 1;
}

void sam_show_header(sam_hdr_t *h) {
	HD_t *hd = &h->hd;
	fprintf(stderr, "@HD:\tVN:%s\tSO:%s\n", hd->VN, hd->SO);
	int i;
	for(i = 0; i < h->sqv.n; ++i) {
		SQ_t *sq = &h->sqv.a[i];
		fprintf(stderr,"@SQ:\tSN:%s\tLN:%d\n", sq->SN, sq->LN);
	}
	PG_t *pg = &h->pg;
	fprintf(stderr,"@PG:\tID:%s\tPN:%s\tVN:%s\tCL:%s\n", pg->ID, pg->PN, pg->VN, pg->CL);
}

static inline int _2char_equal(char a, char b, const char *c) {
	return (a == c[0] && b == c[1]);
}

void sam_header(const char *line, int len, sam_hdr_t *h) {
	int i, cnt_d = 0, l = 0;
	const char delimiter = '\t';
	char type[5], buf[1024];
	SQ_t sq; memset(&sq, 0, sizeof(sq));
	for(i = 0; i <= len; ++i) {
		if(i == len || line[i] == delimiter) {
			buf[l] = '\0'; l = 0;
			++cnt_d;
			if(cnt_d == 1) {
				strcpy(type, buf);
			} else {
				if(!strcmp(type, "HD")) {
					if(_2char_equal(buf[0], buf[1], "VN")) {
						h->hd.VN = strdup(buf + 3);
					} else if(_2char_equal(buf[0], buf[1], "SO")) {
						h->hd.SO = strdup(buf + 3);
					}
				} else if(!strcmp(type, "SQ")) {
					if(_2char_equal(buf[0], buf[1], "SN")) {
						sq.SN = strdup(buf + 3);
					} else if(_2char_equal(buf[0], buf[1], "LN")) {
						sq.LN = atoi(buf + 3);
					}
				} else if(!strcmp(type, "PG")) {
					if(_2char_equal(buf[0], buf[1], "ID")) {
						h->pg.ID = strdup(buf + 3);
					} else if(_2char_equal(buf[0], buf[1], "PN")) {
						h->pg.PN = strdup(buf + 3);
					} else if(_2char_equal(buf[0], buf[1], "VN")) {
						h->pg.VN = strdup(buf + 3);
					} else if(_2char_equal(buf[0], buf[1], "CL")) {
						h->pg.CL = strdup(buf + 3);
					}
				}
			}
		} else {
			buf[l++] = line[i];
		}
	}
	if(!strcmp(type, "SQ")) {
		kv_push(SQ_t, h->sqv, sq);
	}
}

static inline void get_opt(const char *s, sam_opt_t *opt) {
	if(_2char_equal(s[0], s[1], "NM")) {
		opt->nm = atoi(s + 5); // NM:i:
	} else if(_2char_equal(s[0], s[1], "AS")) {
		opt->as = atoi(s + 5); // AS:i:
	}
}

void sam_record1(const char *line, int len, sam_core1_t *r) {
	memset(r, 0, sizeof(*r));
	int i, cnt_d = 0;
	r->data = strdup(line);
	for(i = 0; i < len; ++i) {
		if(r->data[i] == '\t') {
			r->data[i] = '\0';
		}
	}
	for(i = 0; i < len; ++i) {
		if(i == 0 || r->data[i-1] == '\0') {
			if(cnt_d == 0) {
				r->qname = r->data + i;
			} else if(cnt_d == 1) {
				r->flag = atoi(r->data + i);
			} else if(cnt_d == 2) {
				r->rname = r->data + i;
			} else if(cnt_d == 3) {
				r->pos = atoi(r->data + i);
			} else if(cnt_d == 4) {
				r->mapq = atoi(r->data + i);
			} else if(cnt_d == 5) {
				r->cigar = r->data + i;
			} else if(cnt_d == 6) {
				r->rnext = r->data + i;
			} else if(cnt_d == 7) {
				r->pnext = atoi(r->data + i);
			} else if(cnt_d == 8) {
				r->tlen = atoi(r->data + i);
			} else if(cnt_d == 9) {
				r->seq = r->data + i;
			} else if(cnt_d == 10) {
				r->qual = r->data + i;
			} else {
				get_opt(r->data + i, &r->opt);
			}
		}
		if(r->data[i] == '\0') {
			++cnt_d;
		}
	}
}

static long ref_len(const sam_hdr_t *h) {
	long res = 0;
	int i;
	const SQ_v *sqs = &h->sqv;
	for(i = 0; i < sqs->n; ++i) {
		res += sqs->a[i].LN;
	}
	return res;
}

sam_info_t sam_all_records(FILE *f) {
	sam_info_t info;
	sam_hdr_t *h = &info.header; memset(h, 0, sizeof(*h));
	sam_core1_v *s1 = &info.s1, *s2 = &info.s2, *s0 = &info.s0;
	kv_init(*s1); kv_init(*s2); kv_init(*s0);
	char line[65536]; // for next generation sequence, 64K buffer is safe and enough.
	sam_core1_t core;
	size_t r_bytes = 0;
	int cnt_h = 0;
	double r_time = realtime();
	while(fgets(line, sizeof(line), f) != NULL) {
		int len = (int)strlen(line);
		r_bytes += len;
		if(line[len-1] == '\n') { line[--len] = '\0'; }
		if(line[0] == '@') {
			++cnt_h;
			sam_header(line + 1, len - 1, h);
		} else {
			sam_record1(line, len, &core);
			if(is_read1(core.flag) && !is_read2(core.flag)) {
				kv_push(sam_core1_t, *s1, core);
			} else if(is_read2(core.flag) && !is_read1(core.flag)) {
				kv_push(sam_core1_t, *s2, core);
			} else {
				kv_push(sam_core1_t, *s0, core);
			}
		}
	}
	if(info.s1.n != 0 && info.s2.n != 0) {
		assert(info.s0.n == 0); // record belongs to read1 or read2 in paired-end mode.
		info.mode_pe = 1;
	} else if(info.s0.n != 0){
		assert(info.s1.n == 0 && info.s2.n == 0); // single-end mode.
		info.mode_pe = 0;
	} else {
		info.mode_pe = -1;
	}
	assert(info.mode_pe != -1);
	info.ref_len = ref_len(h);
	fprintf(stderr, "[%s] reading time cost %.3f sec, %d header, %ld bytes\n", __func__, realtime()-r_time, cnt_h, r_bytes);
	fprintf(stderr, "    reference length: %ld\n", info.ref_len);
	if(info.mode_pe == 1) {
		fprintf(stderr, "    %ld READ1 records, %ld READ2 records\n", s1->n, s2->n);
	} else {
		fprintf(stderr, "    %ld single-end records\n", s0->n);
	}
	return info;
}

static void coverage1(const sam_core1_v *s, long ref_l) {
	int i;
	long bases = 0;
	for(i = 0; i < s->n; ++i) {
		const sam_core1_t *p = &s->a[i];
		if(!unmap(p->flag)) {
			bases += strlen(p->seq);
		}
	}
	fprintf(stderr, "[%s] %ld bases, %ld reference, ~%.3f coverage\n", __func__, bases, ref_l, 1.0 * bases / ref_l);
}

static void coverage2(const sam_core1_v *s1, const sam_core1_v *s2, long ref_l) {
	int i;
	long bases = 0;
	for(i = 0; i < s1->n; ++i) {
		const sam_core1_t *p = &s1->a[i];
		if(!unmap(p->flag)) {
			bases += strlen(p->seq);
		}
	}
	for(i = 0; i < s2->n; ++i) {
		const sam_core1_t *p = &s2->a[i];
		if(!unmap(p->flag)) {
			bases += strlen(p->seq);
		}
	}
	fprintf(stderr, "[%s] %ld bases, %ld reference, ~%.3f coverage\n", __func__, bases, ref_l, 1.0 * bases / ref_l);
}

static int get_dis(const sam_core1_t *s1, const sam_core1_t *s2) {
	if(strcmp(s1->rname, s2->rname) != 0) {
		return -1;
	} else {
		assert(s2->pos >= s1->pos);
		return s2->pos - s1->pos;
	}
}

static int get_mis(const char *s1, const char *s2, int d) {
	int i, ret = 0, len1 = strlen(s1), len2 = strlen(s2);
	for(i = 0; i < len2 && i+d < len1; ++i) {
		if(s1[i + d] != s2[i]) {
			++ret;
		}
	}
	return ret;
}

void sam_destroy(sam_info_t *info) {
	sam_hdr_t *h = &info->header;
	int i;
	for(i = 0; i < h->sqv.n; ++i) {
		free(h->sqv.a[i].SN);
	}
	free(h->sqv.a);
	free(h->hd.SO); free(h->hd.VN);
	free(h->pg.ID); free(h->pg.PN); free(h->pg.VN); free(h->pg.CL);
	sam_core1_v *s0 = &info->s0, *s1 = &info->s1, *s2 = &info->s2;
	if(s0->n > 0) {
		for(i = 0; i < s0->n; ++i) {
			free(s0->a[i].data);
		}
		free(s0->a);
	}
	if(s1->n > 0) {
		for(i = 0; i < s1->n; ++i) {
			free(s1->a[i].data);
		}
		free(s1->a);
	}
	if(s2->n > 0) {
		for(i = 0; i < s2->n; ++i) {
			free(s2->a[i].data);
		}
		free(s2->a);
	}
}

#define dis_lt(a, b) (strcmp((a).rname,(b).rname) != 0 ?strcmp((a).rname, (b).rname) < 0 :(a).pos < (b).pos)
KSORT_INIT(dis, sam_core1_t, dis_lt)
static void distance2(const sam_core1_v *s1, const sam_core1_v *s2) {
	sam_core1_v s; kv_init(s);
	s.m = s1->n + s2->n;
	s.a = malloc(s.m * sizeof(sam_core1_t));
	int i;
	for(i = 0; i < s1->n; ++i) {
		const sam_core1_t *p = &s1->a[i];
		if(unmap(p->flag) || sec_ali(p->flag) || sup_ali(p->flag)) {
			continue;
		}
		kv_push(sam_core1_t, s, *p);
	}
	for(i = 0; i < s2->n; ++i) {
		const sam_core1_t *p = &s2->a[i];
		if(unmap(p->flag) || sec_ali(p->flag) || sup_ali(p->flag)) {
			continue;
		}
		kv_push(sam_core1_t, s, *p);
	}
	ks_introsort(dis, s.n, s.a);
	long all_d = 0, cnt = 0, mis = 0;
	for(i = 1; i < s.n; ++i) {
		int d = get_dis(&s.a[i-1], &s.a[i]);
		if(d != -1) {
			all_d += d;
			++cnt;
			mis += get_mis(s.a[i-1].seq, s.a[i].seq, d);
		}
	}
	fprintf(stderr, "[%s] Average distance between contiguous reads: %.3f\n", __func__, 1.0 * all_d / cnt);
	fprintf(stderr, "    average mismatches between contiguous reads: %.3f\n", 1.0 * mis / cnt);
	free(s.a);
}

static void distance1(const sam_core1_v *s1) {
	sam_core1_v s; kv_init(s);
	s.m = s1->n;
	s.a = malloc(s.m * sizeof(sam_core1_t));
	int i;
	for(i = 0; i < s1->n; ++i) {
		const sam_core1_t *p = &s1->a[i];
		if(unmap(p->flag) || sec_ali(p->flag) || sup_ali(p->flag)) {
			continue;
		}
		kv_push(sam_core1_t, s, *p);
	}
	ks_introsort(dis, s.n, s.a);
	long all_d = 0, cnt = 0, mis = 0;
	for(i = 1; i < s.n; ++i) {
		int d = get_dis(&s.a[i-1], &s.a[i]);
		if(d != -1) {
			all_d += d;
			++cnt;
			mis += get_mis(s.a[i-1].seq, s.a[i].seq, d);
		}
	}
	fprintf(stderr, "[%s] Average distance between contiguous reads: %.3f\n", __func__, 1.0 * all_d / cnt);
	fprintf(stderr, "    average mismatches between contiguous reads: %.3f\n", 1.0 * mis / cnt);
	free(s.a);
}

static void statistics1(const sam_core1_v *v) {
	int i;
	int unmap_n = 0;
	int sec_n = 0, pri_n = 0, sup_n = 0;
	for(i = 0; i < v->n; ++i) {
		const sam_core1_t *r = &v->a[i];
		if(unmap(r->flag)) ++unmap_n;
		if(sec_ali(r->flag)) ++sec_n;
		if(sup_ali(r->flag)) ++sup_n;
		if(!sec_ali(r->flag) && !sup_ali(r->flag)) ++pri_n;
	}
	fprintf(stderr, "[%s] Mapping statistics for single-end reads\n", __func__);
	fprintf(stderr, "    # records          %ld\n", v->n);
	fprintf(stderr, "    # unmap            %d ~ %.3f %%\n", unmap_n, 100.0 * unmap_n / v->n);
	fprintf(stderr, "    # secondary        %d ~ %.3f %%\n", sec_n, 100.0 * sec_n / v->n);
	fprintf(stderr, "    # supplementary    %d ~ %.3f %%\n", sup_n, 100.0 * sup_n / v->n);
	fprintf(stderr, "    # primary          %d ~ %.3f %%\n", pri_n, 100.0 * pri_n / v->n);
}

static void statistics2(const sam_core1_v *v1, const sam_core1_v *v2) {
	int i;
	int unmap1_n = 0,  both_map1 = 0;
	int sec1_n = 0, pri1_n = 0, sup1_n = 0;
	for(i = 0; i < v1->n; ++i) {
		const sam_core1_t *r = &v1->a[i];
		if(unmap(r->flag)) ++unmap1_n;
		if(both_ali(r->flag)) ++both_map1;
		if(sec_ali(r->flag)) ++sec1_n;
		if(sup_ali(r->flag)) ++sup1_n;
		if(!sec_ali(r->flag) && !sup_ali(r->flag)) ++pri1_n;
	}

	int unmap2_n = 0, both_map2 = 0;
	int sec2_n = 0, pri2_n = 0, sup2_n = 0;
	for(i = 0; i < v2->n; ++i) {
		const sam_core1_t *r = &v2->a[i];
		if(unmap(r->flag)) ++unmap2_n;
		if(both_ali(r->flag)) ++both_map2;
		if(sec_ali(r->flag)) ++sec2_n;
		if(sup_ali(r->flag)) ++sup2_n;
		if(!sec_ali(r->flag) && !sup_ali(r->flag)) ++pri2_n;
	}
	assert(both_map1 == both_map2);

	fprintf(stderr, "[%s] Mapping statistics for single-end reads\n", __func__);
	fprintf(stderr, "    READ1\n");
	fprintf(stderr, "        # records          %ld\n", v1->n);
	fprintf(stderr, "        # unmap            %d ~ %.3f %%\n", unmap1_n, 100.0 * unmap1_n / v1->n);
	fprintf(stderr, "        # secondary        %d ~ %.3f %%\n", sec1_n, 100.0 * sec1_n / v1->n);
	fprintf(stderr, "        # supplementary    %d ~ %.3f %%\n", sup1_n, 100.0 * sup1_n / v1->n);
	fprintf(stderr, "        # primary          %d ~ %.3f %%\n", pri1_n, 100.0 * pri1_n / v1->n);
	fprintf(stderr, "    READ2\n");
	fprintf(stderr, "        # records          %ld\n", v2->n);
	fprintf(stderr, "        # unmap            %d ~ %.3f %%\n", unmap2_n, 100.0 * unmap2_n / v2->n);
	fprintf(stderr, "        # secondary        %d ~ %.3f %%\n", sec2_n, 100.0 * sec2_n / v2->n);
	fprintf(stderr, "        # supplementary    %d ~ %.3f %%\n", sup2_n, 100.0 * sup2_n / v2->n);
	fprintf(stderr, "        # primary          %d ~ %.3f %%\n", pri2_n, 100.0 * pri2_n / v2->n);
}

static void wgsim_evaluate1(const sam_core1_v *v) {
	const int BWA_MEM_QF_COEF = 3; // quality filter lowerbound
	const int POS_DIFF = 5;
	int i, j, reads_n = 0;
	int wrong[300], mapped[300];
	memset(wrong, 0, sizeof(wrong));
	memset(mapped, 0, sizeof(mapped));
	for(i = 0; i < v->n; i++) {
		const sam_core1_t *p = &v->a[i];
		if(unmap(p->flag)) {
			reads_n++;
			continue;
		}
		if(sup_ali(p->flag) || sec_ali(p->flag)) continue;
		reads_n++;
		if(p->mapq <= BWA_MEM_QF_COEF) continue;
		// MAPQ: [0, 255]
		mapped[p->mapq]++;
		int qn_len = strlen(p->qname);
		int _cnt = 0;
		kstring_t tmp; memset(&tmp, 0, sizeof(tmp));
		int sim_pos1 = -1, sim_pos2 = -1;
		for(j = 0; j < qn_len; j++) {
			if(p->qname[j] == '_') {
				_cnt++;
				if(_cnt == 2) {
					sim_pos1 = atoi(tmp.s);
				} else if(_cnt == 3) {
					sim_pos2 = atoi(tmp.s);
					break;
				}
				tmp.l = 0;
			} else {
				kputc(p->qname[j], &tmp);
			}
		}
		free(tmp.s);
		sim_pos2 = sim_pos2 - strlen(p->seq) + 1;
		if(is_rc(p->flag)) {
			if(abs(p->pos - sim_pos2) > POS_DIFF) {
				wrong[p->mapq]++;
			}
		} else {
			if(abs(p->pos - sim_pos1) > POS_DIFF) {
				wrong[p->mapq]++;
			}
		}
	}
	int cnt = 0;
	for(i = 254; i > BWA_MEM_QF_COEF; i--) {
		mapped[i] += mapped[i+1];
		wrong[i] += wrong[i+1];
	}
	printf("Reads: %d , Wrong: %d , Mapped: %d\n", reads_n, wrong[BWA_MEM_QF_COEF+1], mapped[BWA_MEM_QF_COEF+1]);
	printf("X ---- Mapped / reads_n\n");
	for(i = BWA_MEM_QF_COEF + 1; i < 255; i++) {
		if(mapped[i] == 0) break;
		printf("%.9f, ", 100.0 * mapped[i] / reads_n);
		cnt++;
		if(cnt % 8 == 0) printf("\n");
	}
	printf("\n");

	cnt = 0;
	printf("Y ---- Wrong / mapped\n");
	for(i = BWA_MEM_QF_COEF + 1; i < 255; i++) {
		if(mapped[i] == 0) break;
		printf("%.9f, ", 1.0 * wrong[i] / mapped[i]);
		cnt++;
//		printf(" >= MAPQ: %d , M: %d , W: %d    ", i, mapped[i], wrong[i]);
		if(cnt % 8 == 0) printf("\n");
	}
	printf("\n");
}

int samop_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	int c, lo_index;
	char *in_fn = NULL;
	const char *short_opts = "i:r:";
	const struct option long_opts[] = {
		// {name, has_arg, flag, val}
		{"cov", 0, NULL, 0},
		{"dis", 0, NULL, 0},
		{"stat", 0, NULL, 0},
		{"cs", 0, NULL, 0},
		{"wgsim-eval", 0, NULL, 0},
		{NULL, 0, NULL, 0}
	};

	int show_r = -1;
	int cov = 0, dis = 0, cs = 0, stat = 0, wgsim_eval = 0;
	while((c=getopt_long(argc, argv, short_opts, long_opts, &lo_index)) >= 0) {
		switch (c) {
			case 0:
				if(!strcmp(long_opts[lo_index].name, "cov")) {
					cov = 1;
				}
				else if(!strcmp(long_opts[lo_index].name, "cs")) {
					cs = 1;
				}
				else if(!strcmp(long_opts[lo_index].name, "dis")) {
					dis = 1;
				}
				else if(!strcmp(long_opts[lo_index].name, "stat")) {
					stat = 1;
				}
				else if(!strcmp(long_opts[lo_index].name, "wgsim-eval")) {
					wgsim_eval = 1;
				}
				break;
			case 'i': in_fn = strdup(optarg); break;
			case 'r': show_r = atoi(optarg); break;
			case '?': return usage();
			default : return usage();
		}
	}
	if(in_fn == NULL) { fprintf(stderr, "input file can't be NULL\n"); return usage(); }

	FILE *fi = fopen(in_fn, "r");
	sam_info_t info = sam_all_records(fi); // todo: load all of records costs too much memory. set a buffer instead

	if(info.mode_pe == 1) {
		if(stat) statistics2(&info.s1, &info.s2);
		if(show_r != -1) { fprintf(stderr, "[%s] PE record [%d]\n", __func__, show_r); }
		if(show_r >= 0 && show_r < info.s1.n) { sam_show_record1(&info.s1.a[show_r]); }
		if(show_r >= 0 && show_r < info.s2.n) { sam_show_record1(&info.s2.a[show_r]); }
		if(cov) coverage2(&info.s1, &info.s2, info.ref_len);
		if(dis) distance2(&info.s1, &info.s2);
	} else {
		if(stat) statistics1(&info.s0);
		if(show_r != -1) { fprintf(stderr, "[%s] SE record [%d]\n", __func__, show_r); }
		if(show_r >= 0 && show_r < info.s0.n) { sam_show_record1(&info.s0.a[show_r]); }
		if(cov) coverage1(&info.s0, info.ref_len);
		if(dis) distance1(&info.s0);
		if(wgsim_eval) wgsim_evaluate1(&info.s0);
	}
	sam_destroy(&info);
	return 0;
}

