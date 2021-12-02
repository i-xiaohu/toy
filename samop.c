//
// Created by ixiaohu on 2019/9/19.
//

#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "samop.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    samop handle SAM file\n");
	fprintf(stderr, "Usage:      samop <in.sam.gz>\n");
	fprintf(stderr, "\n");
	return 1;
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

void sam_record1(char *line, int len, sam_core1_t *r) {
	memset(r, 0, sizeof(*r));
	int i, cnt_d = 0;
	r->data = line;
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

static long get_ref_len(const sam_hdr_t *h) {
	long res = 0;
	int i;
	const SQ_v *sqs = &h->sqv;
	for(i = 0; i < sqs->n; ++i) {
		res += sqs->a[i].LN;
	}
	return res;
}

static void sam_header_free(sam_hdr_t *h) {
	int i;
	for(i = 0; i < h->sqv.n; ++i) {
		free(h->sqv.a[i].SN);
	}
	free(h->sqv.a);
	free(h->hd.SO); free(h->hd.VN);
	free(h->pg.ID); free(h->pg.PN); free(h->pg.VN); free(h->pg.CL);
}

static void analysis(gzFile f) {
	sam_hdr_t h; memset(&h, 0, sizeof(h));
	char line[65536]; // for next generation sequence, 64K buffer is safe and enough.
	sam_core1_t core;
	long cnt = 0;
	long unmap_n = 0;
	long sec_n = 0, pri_n = 0, sup_n = 0;
	long as = 0, nm = 0;
	while(gzgets(f, line, sizeof(line)) != NULL) {
		int len = (int)strlen(line);
		if(line[len-1] == '\n') { line[--len] = '\0'; }
		if(line[0] == '@') {
			sam_header(line + 1, len - 1, &h); // skip the char '@'
		} else {
			sam_record1(line, len, &core);
			const sam_core1_t *r = &core;
			cnt++;
			if(unmap(r->flag)) {
				unmap_n++;
				continue;
			}
			if(sec_ali(r->flag)) sec_n++;
			if(sup_ali(r->flag)) sup_n++;
			if(!sec_ali(r->flag) && !sup_ali(r->flag)) {
				pri_n++;
				as += r->opt.as;
				nm += r->opt.nm;
			}
		}
	}

	long ref_len = get_ref_len(&h);
	sam_header_free(&h);
	fprintf(stderr, "Reference_length: %ld\n", ref_len);
	fprintf(stderr, "Alignments:       %ld\n", cnt);
	fprintf(stderr, "Unmap:            %ld\n", unmap_n);
	fprintf(stderr, "Primary:          %ld\n", pri_n);
	fprintf(stderr, "Supplementary:    %ld\n", sup_n);
	fprintf(stderr, "Secondary:        %ld\n", sec_n);
	fprintf(stderr, "Alignment_Score:  %.2f\n", 1.0 * as / pri_n);
	fprintf(stderr, "Edit_distance:    %.2f\n", 1.0 * nm / pri_n);
	fprintf(stderr, "\n");
}

int samop_main(int argc, char *argv[]) {
	if(argc != 2) {
		return usage();
	}

	gzFile fi = gzopen(argv[1], "r");
	if (fi == NULL) {
		fprintf(stderr, "Open %s failed\n", argv[1]);
		return 1;
	}
	analysis(fi);
	gzclose(fi);
	return 0;
}

