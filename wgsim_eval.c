//
// Created by ixiaohu on 2021/1/15.
//

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include "samop.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    weval evaluate the SAM file of wgsim reads\n");
	fprintf(stderr, "Usage:      weval [options] <in.sam>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "            -d [INT] The maximal distance allowed [10bp].\n");
	fprintf(stderr, "Note: only primary alignments are considered.\n");
	fprintf(stderr, "\n");
	return 1;
}

typedef struct {
	int dis;
	int lines_n;
	FILE *f;
} opt_t;

static opt_t* opt_init() {
	opt_t* ret = malloc(sizeof(opt_t));
	ret->dis = 10;
	ret->lines_n = 1000*1000; // 1M SAM lines.
	ret->f = NULL;
	return ret;
}

static void opt_destroy(opt_t *opt) {
	fclose(opt->f);
	free(opt);
}

typedef struct {
	opt_t *opt;
	int correct[255];
	int align[255];
	int miss_n; // #unmapped reads
} ktp_aux_t;

static sam_core1_v* tp_input(FILE *f, int lines_n) {
	char line[65536]; // for next generation sequence, 64K buffer is safe and enough.
	int cnt = 0;
	sam_core1_v *ret = calloc(1, sizeof(sam_core1_v));
	while (fgets(line, sizeof(line), f) != NULL) {
		if(line[0] == '@') continue;
		cnt++;
		int len = strlen(line);
		if(line[len-1] == '\n') { line[--len] = '\0'; }
		sam_core1_t r;
		sam_record1(line, len, &r);
		kv_push(sam_core1_t, *ret, r);
		if(cnt >= lines_n) break;
	}
	return ret;
}

static void tp_process(sam_core1_v *data, ktp_aux_t *aux) {
	int i, j;
	opt_t *opt = aux->opt;
	for(i = 0; i < data->n; i++) {
		const sam_core1_t *p = &data->a[i];
		if(unmap(p->flag)) {
			aux->miss_n++;
			continue;
		}
		if(sec_ali(p->flag) || sup_ali(p->flag)) continue; // Only primary alignments considered.
		aux->align[p->mapq]++;
		// Extract the simulated position from query name in format 'rname_read1Pos_read2Pos'
		char rname[64]; memset(rname, 0, sizeof(rname));
		for(j=0; p->qname[j] && p->qname[j]!='_'; j++) {
			rname[j] = p->qname[j];
		}
		if(!p->qname[j]) {
			fprintf(stderr, "Reference-name absent in wgsim read name.\n");
			abort();
		}
		int sim_pos1 = 0, sim_pos2 = 0;
		for(j=j+1; p->qname[j] && p->qname[j]!='_'; j++) {
			sim_pos1 *= 10; sim_pos1 += p->qname[j]-'0';
		}
		if(!p->qname[j]) {
			fprintf(stderr, "read1Pos absent in wgsim read name.\n");
			abort();
		}
		for(j=j+1; p->qname[j] && p->qname[j]!='_'; j++) {
			sim_pos2 *= 10; sim_pos2 += p->qname[j]-'0';
		}
		if(!p->qname[j]) {
			fprintf(stderr, "read2Pos absent in wgsim read name.\n");
			abort();
		}
		int qlen = strlen(p->seq);
		sim_pos2 = sim_pos2 - qlen + 1; // Correct the read2Pos on the forward strand.
		if(!is_rc(p->flag) && abs(p->pos-sim_pos1) <= opt->dis) aux->correct[p->mapq]++;
		if( is_rc(p->flag) && abs(p->pos-sim_pos2) <= opt->dis) aux->correct[p->mapq]++;
	}
}

static void tp_free(sam_core1_v *data) {
	int i;
	for(i = 0; i < data->n; i++) {
		free(data->a[i].data);
	}
	free(data->a);
	free(data);
}

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

static void* tp_main(void *_aux, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)_aux;
	sam_core1_v *data = (sam_core1_v*)_data;
	opt_t *opt = aux->opt;
	if(step == 0) {
		data = tp_input(opt->f, opt->lines_n);
		if(data->n == 0) { free(data); return NULL; }
		fprintf(stderr, "Input %ld records\n", data->n);
		return data;
	}
	else if(step == 1) {
		tp_process(data, aux);
		return data;
	}
	else if(step == 2){
		tp_free(data);
		return NULL;
	}
	return NULL;
}

int weval_main(int argc, char** argv) {
	if(argc == 1) return usage();
	int c;
	opt_t *opt = opt_init();
	while ((c = getopt(argc, argv, "d:")) >= 0) {
		if (c == 'd') opt->dis = atoi(optarg);
	}
	if(optind != argc-1) return usage();

	opt->f = fopen(argv[optind], "r");
	if(opt->f == NULL) {
		fprintf(stderr, "Open %s filed.\n", argv[optind]);
		return 1;
	}

	ktp_aux_t *aux = calloc(1, sizeof(ktp_aux_t));
	aux->opt = opt;
	kt_pipeline(2, tp_main, aux, 3);
	fprintf(stdout, "%d\n", aux->miss_n);

	int i;
	for (i = 0; i < 255; i++) {
		if(aux->align[i] == 0) continue;
		fprintf(stdout, "%d %d %d\n", i, aux->align[i], aux->correct[i]);
	}

	opt_destroy(opt);
	return 0;
}