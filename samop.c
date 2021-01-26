//
// Created by ixiaohu on 2019/9/19.
//

#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "samop.h"
#include "utils.h"
#include "ksort.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    samop handle SAM file\n");
	fprintf(stderr, "Usage:      samop [options] <in.sam>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "    -r INT  Show a (pair of) SAM record\n");
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
	fprintf(stderr, "\n");
	fprintf(stderr, "Mapping statistics for single-end reads\n");
	fprintf(stderr, "    #records          %ld\n", v->n);
	fprintf(stderr, "    #unmap            %d\n", unmap_n);
	fprintf(stderr, "    #primary          %d\n", pri_n);
	fprintf(stderr, "    #supplementary    %d\n", sup_n);
	fprintf(stderr, "    #secondary        %d\n", sec_n);
	fprintf(stderr, "\n");
}

int samop_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}

	int c, show_r = -1;
	while((c=getopt(argc, argv, "r:")) >= 0) {
		if(c == 'r') {
			show_r = atoi(optarg);
		}
		else return usage();
	}

	if(optind != argc-1) {
		return usage();
	}

	FILE *fi = fopen(argv[optind], "r");
	sam_info_t info = sam_all_records(fi); // todo: load all of records costs too much memory. set a buffer instead

	if(info.mode_pe == 1) {
		fprintf(stderr, "Paired-end alignments are not supportive.\n");
	} else {
		statistics1(&info.s0);
		if(show_r != -1) fprintf(stderr, "Show SE record [%d]\n", show_r);
		if(show_r >= 0 && show_r < info.s0.n) sam_show_record1(&info.s0.a[show_r]);
		else fprintf(stderr, "%d is out of the boundary.\n", show_r);
	}
	sam_destroy(&info);
	return 0;
}

