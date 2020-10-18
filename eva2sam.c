//
// Created by 63175 on 2019/11/1.
//

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "samop.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:    eva2sam   evaluate the similarity of two alignment files\n");
	fprintf(stderr, "Usage:      eva2sam   [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "            -z STR        MEZ file name\n");
	fprintf(stderr, "            -m STR        MEM file name\n");
	fprintf(stderr, "            -k INT        buffer size for each SAM records\n");
	fprintf(stderr, "\n");
	return 1;
}

typedef struct {
	int z0m0, z1m0, z0m1, z1m1; // the count of mapped/unmapped reads
	int z1m0_has_N;
	int h_z1m1, h_z0m1, h_z1m0; // supplementary or secondary alignment
	int nm_zlm, nm_zgm, nm_zem; // NM compare statistics
	int nm3_zlm, nm3_zgm, nm3_zem; // only keep records MAPQ greater than 3 (BWA-MEM paper does it)
	int pos3_zem, pos3_znem;
	int as_zlm, as_zem, as_zgm; // alignment score comparison
	int pri_as_zlm, pri_as_zem, pri_as_zgm; // alignment score comparison only for primary alignments
} stat_t;

/* only keep one SAM record  */
static void do_check(FILE *fz, FILE *fm) {
	sam_hdr_t mezh; memset(&mezh, 0, sizeof(mezh));
	sam_hdr_t memh; memset(&memh, 0, sizeof(memh));
	char lm_line[65536], lz_line[65536];

	/* MEM header */
	while(fgets(lm_line, sizeof(lm_line), fm) != NULL) {
		int len = (int)strlen(lm_line);
		if(lm_line[len-1] == '\n') { lm_line[--len] = '\0'; }
		if(lm_line[0] == '@') {
			sam_header(lm_line + 1, len - 1, &memh);
		} else {
			break;
		}
	}
	if((!memh.hd.SO) || (strcmp(memh.hd.SO, "queryname") != 0)) {
		fprintf(stderr, "Warning: MEM SAM isn't sorted by query name\n");
	}

	/* MEZ header */
	while(fgets(lz_line, sizeof(lz_line), fz) != NULL) {
		int len = (int)strlen(lz_line);
		if(lz_line[len-1] == '\n') { lz_line[--len] = '\0'; }
		if(lz_line[0] == '@') {
			sam_header(lz_line + 1, len - 1, &mezh);
		} else {
			break;
		}
	}
	if(!mezh.hd.SO || strcmp(mezh.hd.SO, "queryname") != 0) {
		fprintf(stderr, "Warning: MEZ SAM isn't sorted by query name\n");
	}

	/* check records */
	sam_core1_t zr, next_z; memset(&next_z, 0, sizeof(next_z));
	sam_record1(lz_line, (int)strlen(lz_line), &zr);
	sam_core1_t mr, next_m; memset(&next_m, 0, sizeof(next_m));
	sam_record1(lm_line, (int)strlen(lm_line), &mr);
	stat_t s; memset(&s, 0, sizeof(s));
	int reads_n = 0;
	for( ;; ) {
		++reads_n;
		/* go forward and check one query read each time */
		int z_has_sup = 0, z_has_next = 0, z_as = zr.opt.as;
		while(fgets(lz_line, sizeof(lz_line), fz) != NULL) {
			int lz = (int)strlen(lz_line);
			if(lz_line[lz - 1] == '\n') { lz_line[--lz] = '\0'; }
			sam_record1(lz_line, lz, &next_z);
			if(!strcmp(zr.qname, next_z.qname)) { // has secondary or supplementary alignments
				z_as += next_z.opt.as;
				free(next_z.data);
				z_has_sup = 1;
			} else {
				z_has_next = 1;
				break;
			}
		}

		int m_has_sup = 0, m_has_next = 0, m_as = mr.opt.as;
		while(fgets(lm_line, sizeof(lm_line), fm) != NULL) {
			int lm = (int)strlen(lm_line);
			if(lm_line[lm - 1] == '\n') { lm_line[--lm] = '\0'; }
			sam_record1(lm_line, lm, &next_m);
			if(!strcmp(mr.qname, next_m.qname)) { // has secondary or supplementary alignments
				m_as += next_m.opt.as;
				free(next_m.data);
				m_has_sup = 1;
			} else {
				m_has_next = 1;
				break;
			}
		}
		assert(strcmp(mr.qname, zr.qname) == 0);

		/* compare two SAM records */
		if(unmap(zr.flag) && unmap(mr.flag)) {
			++s.z0m0;
		} else if(!unmap(zr.flag) && unmap(mr.flag)) {
			++s.z1m0;
			int i; // if we aligned with N that BWA-MEM ignored
			for(i = 0; mr.seq[i]; ++i) {
				if(mr.seq[i] == 'N') {
					++s.z1m0_has_N;
					break;
				}
			}
		} else if(unmap(zr.flag) && !unmap(mr.flag)) {
			++s.z0m1;
		} else {
			++s.z1m1; // both aligned
			/* Compare AS for both aligned reads (sum up when chimeric reads) */
			if(z_as < m_as) ++s.as_zlm;
			else if(z_as == m_as) ++s.as_zem;
			else ++s.as_zgm;

			if(m_has_sup && z_has_sup) {
				++s.h_z1m1;
			} else if(!z_has_sup && m_has_sup) {
				++s.h_z0m1;
			} else if(z_has_sup && !m_has_sup) {
				++s.h_z1m0;
			} else {
				// both primary alignments.
				/* When both MAPQ are greater than 3, compare NM and pos */
				if(zr.mapq > 3) {
					if(zr.pos != mr.pos) {
						++s.pos3_znem;
					} else {
						++s.pos3_zem;
					}
					if(zr.opt.nm < mr.opt.nm) {
						++s.nm3_zlm;
					} else if(zr.opt.nm > mr.opt.nm) {
						++s.nm3_zgm;
					} else if(zr.opt.nm == mr.opt.nm) {
						++s.nm3_zem;
					}
				}
				/* Compare NM for both primary alignments. */
				if(zr.opt.nm < mr.opt.nm) {
					++s.nm_zlm;
				} else if(zr.opt.nm > mr.opt.nm) {
					++s.nm_zgm;
				} else if(zr.opt.nm == mr.opt.nm) {
					++s.nm_zem;
				}
				/* Compare AS for both primary alignments. */
				if(zr.opt.as < mr.opt.as) ++s.pri_as_zlm;
				else if(zr.opt.as == mr.opt.as) ++s.pri_as_zem;
				else ++s.pri_as_zgm;
			}
		}
		free(zr.data); free(mr.data);
		if(!z_has_next || !m_has_next) {
			/* get all records corresponding all reads */
			assert(!z_has_next && !m_has_next); // count of reads must be equal, but count of records may be not
			break;
		} else {
			zr = next_z; mr = next_m;
		}
	}
	fprintf(stderr, "MEZ-MEM get reads  %d\n", reads_n);
	fprintf(stderr, "z0m0               %d ~ %.3f %%\n", s.z0m0, 100.0 * s.z0m0 / reads_n);
	fprintf(stderr, "z1m0               %d ~ %.3f %%\n", s.z1m0, 100.0 * s.z1m0 / reads_n);
	fprintf(stderr, "    has_N          %d ~ %.3f %%\n", s.z1m0_has_N, 100.0 * s.z1m0_has_N/ s.z1m0);
	fprintf(stderr, "z0m1               %d ~ %.3f %%\n", s.z0m1, 100.0 * s.z0m1 / reads_n);
	fprintf(stderr, "z1m1               %d ~ %.3f %%\n", s.z1m1, 100.0 * s.z1m1 / reads_n);
	fprintf(stderr, "    Alignment score comparison(sum up when supplementary alignments)\n");
	fprintf(stderr, "        z < m          %d ~ %.6f %%\n", s.as_zlm, 100.0 * s.as_zlm / s.z1m1);
	fprintf(stderr, "        z = m          %d ~ %.6f %%\n", s.as_zem, 100.0 * s.as_zem / s.z1m1);
	fprintf(stderr, "        z > m          %d ~ %.6f %%\n", s.as_zgm, 100.0 * s.as_zgm / s.z1m1);
	fprintf(stderr, "    supplementary  or secondary alignments statistics\n");
	fprintf(stderr, "        z1m0           %d\n", s.h_z1m0);
	fprintf(stderr, "        z0m1           %d\n", s.h_z0m1);
	fprintf(stderr, "        z1m1           %d\n", s.h_z1m1);
	fprintf(stderr, "    Alignment score comparison for primary alignment\n");
	fprintf(stderr, "        z < m          %d ~ %.6f %%\n", s.pri_as_zlm, 100.0 * s.pri_as_zlm / (s.pri_as_zlm + s.pri_as_zem + s.pri_as_zgm));
	fprintf(stderr, "        z = m          %d ~ %.6f %%\n", s.pri_as_zem, 100.0 * s.pri_as_zem / (s.pri_as_zlm + s.pri_as_zem + s.pri_as_zgm));
	fprintf(stderr, "        z > m          %d ~ %.6f %%\n", s.pri_as_zgm, 100.0 * s.pri_as_zgm / (s.pri_as_zlm + s.pri_as_zem + s.pri_as_zgm));
	fprintf(stderr, "    NM statistics for primary alignments\n");
	fprintf(stderr, "        z < m          %d ~ %.6f %%\n", s.nm_zlm, 100.0 * s.nm_zlm / (s.nm_zlm + s.nm_zem + s.nm_zgm));
	fprintf(stderr, "        z = m          %d ~ %.6f %%\n", s.nm_zem, 100.0 * s.nm_zem / (s.nm_zlm + s.nm_zem + s.nm_zgm));
	fprintf(stderr, "        z > m          %d ~ %.6f %%\n", s.nm_zgm, 100.0 * s.nm_zgm / (s.nm_zlm + s.nm_zem + s.nm_zgm));
	fprintf(stderr, "    NM with MAPQ greater than 3 statistics\n");
	fprintf(stderr, "        z < m          %d ~ %.6f %%\n", s.nm3_zlm, 100.0 * s.nm3_zlm / (s.nm3_zlm + s.nm3_zem + s.nm3_zgm));
	fprintf(stderr, "        z = m          %d ~ %.6f %%\n", s.nm3_zem, 100.0 * s.nm3_zem / (s.nm3_zlm + s.nm3_zem + s.nm3_zgm));
	fprintf(stderr, "        z > m          %d ~ %.6f %%\n", s.nm3_zgm, 100.0 * s.nm3_zgm / (s.nm3_zlm + s.nm3_zem + s.nm3_zgm));
	fprintf(stderr, "    RPOS with MAPQ greater than 3 statistics\n");
	fprintf(stderr, "        z = m          %d ~ %.6f %%\n", s.pos3_zem,  100.0 * s.pos3_zem / (s.pos3_znem + s.pos3_zem));
	// todo: if distance difference isn't greater than 20 after clipping, we think they are equal just like BWA-MEM evaluation scheme
	fprintf(stderr, "        z != m         %d ~ %.6f %%\n", s.pos3_znem, 100.0 * s.pos3_znem / (s.pos3_znem + s.pos3_zem));
}

int eva2sam_main(int argc, char *argv[]) {
	if(argc == 1) {
		return usage();
	}
	const char *opts = "z:m:";
	int c;
	FILE *fz = NULL, *fm = NULL;
	while((c=getopt(argc, argv, opts)) > -1){
		if(c == 'z') {
			fz = fopen(optarg, "r");
		} else if(c == 'm') {
			fm = fopen(optarg, "r");
		} else {
			return usage();
		}
	}
	if(fz == NULL || fm == NULL) {
		fprintf(stderr, "fail to open MEZ/MEM SAM file\n");
		return 1;
	}
	do_check(fz, fm);
	return 0;
}
