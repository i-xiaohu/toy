//
// Created by ixiaohu on 2021/1/15.
//

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

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
} opt_t;

static opt_t* opt_init() {
	opt_t* ret = malloc(sizeof(opt_t));
	ret->dis = 10;
	return ret;
}

static void analysis(gzFile f, const opt_t *opt) {
	int i, j;
	char line[65536];
	int correct[255]; memset(correct, 0, sizeof(correct));
	int align[255]; memset(align, 0, sizeof(align));
	int miss_n = 0; // #unmapped reads
	sam_core1_t core;
	while(gzgets(f, line, sizeof(line)) != NULL) {
		int len = (int)strlen(line);
		if(line[len-1] == '\n') { line[--len] = '\0'; }
		if(line[0] == '@') {
			continue;
		} else {
			sam_record1(line, len, &core);
			if(is_read1(core.flag) || is_read2(core.flag)) {
				fprintf(stderr, "Paired-end alignments are not supportive.\n");
				abort();
			}
			const sam_core1_t *p = &core;
			if(unmap(p->flag)) {
				miss_n++;
				continue;
			}
			if(sec_ali(p->flag) || sup_ali(p->flag)) continue; // Only primary alignments considered.
			align[p->mapq]++;
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
			if(strcmp(rname, p->rname) != 0) continue;
			if(!is_rc(p->flag) && abs(p->pos-sim_pos1) <= opt->dis) correct[p->mapq]++;
			if( is_rc(p->flag) && abs(p->pos-sim_pos2) <= opt->dis) correct[p->mapq]++;
		}
	}
	fprintf(stdout, "%d\n", miss_n);
	for (i = 0; i < 255; i++) {
		if(align[i] == 0) continue;
		fprintf(stdout, "%d %d %d\n", i, align[i], correct[i]);
	}
}

int weval_main(int argc, char** argv) {
	if(argc == 1) return usage();
	int c;
	opt_t *opt = opt_init();
	while ((c = getopt(argc, argv, "d:")) >= 0) {
		if (c == 'd') opt->dis = atoi(optarg);
	}
	if(optind != argc-1) return usage();

	gzFile f = gzopen(argv[optind], "r");
	if(f == NULL) {
		fprintf(stderr, "Open %s filed.\n", argv[optind]);
		return 1;
	}
	analysis(f, opt);

	gzclose(f);
	free(opt);
	return 0;
}