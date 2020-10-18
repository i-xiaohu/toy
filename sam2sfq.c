//
// Created by 63175 on 2019/10/28.
//

#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include "samop.h"
#include "progress.h"
#include "ksort.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    sam2sfq   generate sorted FASTQ file from Coordinate-sorted SAM file\n");
	fprintf(stderr, "                      consecutive reads share the largest overlap in sorted FASTQ file\n");
	fprintf(stderr, "Usage:      sam2sfq   [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -i STR        SAM file name\n");
	fprintf(stderr, "            -o STR        FASTQ file name\n");
	fprintf(stderr, "            -K INT        pairing buffer size [10000000], only in PE mode\n");
	fprintf(stderr, "\n");
	return 1;
}

static void write_sfq1(sam_core1_v *s, FILE *f) {
	if(!s->n) return ;
	fprintf(stderr, "[%s] Output sorted FASTQ...\n", __func__);
	progress_t bar;
	progress_init(&bar, "", 100, PROGRESS_CHR_STYLE);
	int i, j, reads_n = 0, min_l = INT_MAX, max_l = 0;
	for(i = 0; i < s->n; ++i) {
		sam_core1_t *r = &s->a[i];
		if(sup_ali(r->flag) || sec_ali(r->flag)) { // don't output secondary / supplementary segment.
			continue;
		}
		int has_N = 0;
		for(j = 0; r->seq[j]; ++j) {
			if(r->seq[j] == 'N') {
				has_N = 1;
				break;
			}
		}
		if(unmap(r->flag) && has_N) { // reads quality control: throw out unmapped reads with N(has too many N).
			continue;
		}
		++reads_n;
		int len = (int)strlen(r->seq);
		min_l = (min_l < len) ?min_l :len;
		max_l = (max_l > len) ?max_l :len;
		fputs("@", f); fputs(r->qname, f); fputs("\n", f);
		fputs(r->seq, f); fputs("\n", f);
		fputs("+", f); fputs("\n", f);
		fputs(r->qual, f); fputs("\n", f);
		progress_show(&bar, 1.0*(i+1)/s->n);
	}
	progress_destroy(&bar);
	fprintf(stderr, "\n");
	fprintf(stderr, "    done! get %d reads\n", reads_n);
	fprintf(stderr, "    min-read-length %d\n", min_l);
	fprintf(stderr, "    max-read-length %d\n", max_l);
}
static void keep_primary(sam_core1_v *s) {
	int i, old_n = (int)s->n;
	s->n = 0;
	for(i = 0; i < old_n; ++i) {
		sam_core1_t *r = &s->a[i];
		if(sup_ali(r->flag) || sec_ali(r->flag)) {
			free(r->data);
		} else {
			kv_push(sam_core1_t, *s, *r);
		}
	}
}

#define sam_core1_lt(a, b) (strcmp((a).rname, (b).rname) != 0 ?strcmp((a).rname, (b).rname) < 0 :(a).pos < (b).pos)
KSORT_INIT(sam_core1, sam_core1_t, sam_core1_lt)

static void write_sfq2(sam_core1_v *s1, sam_core1_v *s2, FILE *f, int K) {
	if(!s1->n || !s2->n) return ;
	fprintf(stderr, "Discard secondary and supplementary records...\n");
	keep_primary(s1);
	keep_primary(s2);
	assert(s1->n == s2->n);
	fprintf(stderr, "    READ1, READ2 primary records: %ld\n", s1->n);
	fprintf(stderr, "Output sorted FASTQ...\n");
	progress_t bar;
	progress_init(&bar, "", 100, PROGRESS_CHR_STYLE);
	int i, j, pairs_n = 0, min_l = INT_MAX, max_l = 0, bytes = 0;
	sam_core1_v out_v; kv_init(out_v);
	for(i = 0; i < s1->n; ++i) {
		progress_show(&bar, 1.0 * (i + 1) / s1->n);
		sam_core1_t *r1 = &s1->a[i], *r2 = &s2->a[i];
		assert(strcmp(r1->qname, r2->qname) == 0);
		int len1 = (int)strlen(r1->seq), drop1 = 0, has_N = 0;
		for(j = 0; j < len1; ++j) {
			if(r1->seq[j] == 'N') {
				has_N = 1;
				break;
			}
		}
		if(unmap(r1->flag) && has_N) { drop1 = 1; }
		int len2 = (int)strlen(r2->seq), drop2 = 0; has_N = 0;
		for(j = 0; j < len2; ++j) {
			if(r2->seq[j] == 'N') {
				has_N = 1;
				break;
			}
		}
		if(unmap(r2->flag) && has_N) { drop2 = 1; }
		if(drop1 || drop2) { // reads quality control: throw out unmapped reads with N(has too many N).
			continue;
		}
		++pairs_n;
		kv_push(sam_core1_t, out_v, *r1); kv_push(sam_core1_t, out_v, *r2);
		bytes += len1 + len2;
		if(bytes > K) {
			bytes = 0;
			ks_introsort(sam_core1, out_v.n, out_v.a);
			for(j = 0; j < out_v.n; ++j) {
				sam_core1_t *r = &out_v.a[j];
				int len = (int)strlen(r->seq);
				min_l = (min_l < len) ?min_l :len;
				max_l = (max_l > len) ?max_l :len;
				fputs("@", f); fputs(r->qname, f); if(is_read1(r->flag)) fputs("/1\n", f); else fputs("/2\n", f);
				fputs(r->seq, f); fputs("\n", f);
				fputs("+", f); fputs("\n", f);
				fputs(r->qual, f); fputs("\n", f);
			}
			out_v.n = 0;
		}
	}
	progress_destroy(&bar);
	for(j = 0; j < out_v.n; ++j) {
		sam_core1_t *r = &out_v.a[j];
		int len = (int)strlen(r->seq);
		min_l = (min_l < len) ?min_l :len;
		max_l = (max_l > len) ?max_l :len;
		fputs("@", f); fputs(r->qname, f); if(is_read1(r->flag)) fputs("/1\n", f); else fputs("/2\n", f);
		fputs(r->seq, f); fputs("\n", f);
		fputs("+", f); fputs("\n", f);
		fputs(r->qual, f); fputs("\n", f);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "    done! get %d pairs\n", pairs_n);
	fprintf(stderr, "    min-read-length %d\n", min_l);
	fprintf(stderr, "    max-read-length %d\n", max_l);
}

int sam2sortfq_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	const char *opts = "i:o:K:";
	int c, K = 10000000;
	FILE *fi = NULL, *fo = NULL;
	while((c=getopt(argc, argv, opts)) > -1) {
		if(c == 'i') {
			fi = fopen(optarg, "r");
		} else if(c == 'o') {
			fo = fopen(optarg, "w");
		} else if(c == 'K') {
			K = atoi(optarg);
		} else {
			return usage();
		}
	}
	if(fi == NULL) {
		fprintf(stderr, "fail to open SAM file\n");
		return 1;
	}
	sam_info_t info = sam_all_records(fi);
	if(info.mode_pe == 0) {
		fprintf(stderr, "[%s] is working in SE mode\n", __func__);
		if(!info.header.hd.SO || strcmp(info.header.hd.SO, "coordinate") != 0) {
			fprintf(stderr, "Warning: the input SAM isn't sorted by coordinate.\n");
		}
		write_sfq1(&info.s0, fo);
	} else {
		fprintf(stderr, "[%s] is working in PE mode\n", __func__);
		write_sfq2(&info.s1, &info.s2, fo, K);
	}
	return 0;
}
