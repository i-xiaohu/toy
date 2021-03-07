//
// Created by ixiaohu on 2021/1/24.
//

#include <stdio.h>

#include "hfastq.h"
#include "kvec.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: view-ref <ref.fa>\n");
	fprintf(stderr, "Usage: View reference online, input \"RName coordinate[1-based index] length\",");
	fprintf(stderr, " stop viewing with \"exit\".\n");
	fprintf(stderr, "\n");
	return 1;
}

typedef kvec_t(bseq1_t) bseq1_v;
static bseq1_v load_all(const open_fastq_t *of) {
	bseq1_v seqs; kv_init(seqs);
	do {
		bseq1_t x = hfastq_fetch1(of);
		kv_push(bseq1_t, seqs, x);
		if(x.name == NULL) break;
	} while(1);
	fprintf(stderr, "[%s] Loading done\n", __func__ );
	return seqs;
}

static void destroy_seqs(const bseq1_v *v) {
	int i;
	for (i = 0; i < v->n; i++) {
		const bseq1_t *p = &v->a[i];
		free(p->name); free(p->comment);
		free(p->seq); free(p->qual);
	}
	free(v->a);
}

int view_ref_main(int argc, char **argv) {
	if(argc == 1) return usage();
	if(argc != 2) return usage();
	open_fastq_t *of = NULL;
	of = hfastq_open(argv[1]);
	if(of == NULL) return 1;

	bseq1_v seqs = load_all(of);

	while (1) {
		fprintf(stderr, "Input \"rname coordinate length\" to view\n");

		char buf[1024];
		fgets(buf, sizeof(buf), stdin);
		if(strcmp(buf, "exit\n") == 0) break;

		char rname[64];
		int cor, len;
		sscanf(buf, "%s %d %d\n", rname, &cor, &len);
		cor--; // 0-based index

		int i, j;
		for (i = 0; i < seqs.n; i++) {
			const bseq1_t *p = &seqs.a[i];
			if(strcmp(rname, p->name) == 0) {
				printf("<%s, %s> Length: %d\n", p->name, p->comment, p->l_seq);
				if(cor < 0 || cor + len > p->l_seq) {
					printf("Index out of the boundary.\n");
					break;
				}
				for (j = 0; j < len; j++) {
					printf("%c", p->seq[cor + j]);
				}
				printf("\n");
				break;
			}
		}
		printf("\n");
	}

	destroy_seqs(&seqs);
	hfastq_close(of);
	return 0;
}