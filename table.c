//
// Created by ixiaohu on 2019/12/1.
//

#include <malloc.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>

#include "table.h"
#include "ksort.h"

tab_t *tab_init(int n, int m, char *header, char *info) {
	tab_t *t = calloc(1, sizeof(tab_t));
	if(header) t->header = strdup(header);
	if(info) t->info = strdup(info);
	t->n = n; t->m = m;
	int i;
	t->cells = malloc(t->n * sizeof(char**));
	for(i = 0; i < n; ++i) {
		t->cells[i] = calloc(t->m, sizeof(char*));
	}
	return t;
}


void tab_fill_cell1(tab_t *t, int x, int y, const char *fmt, ...) {
	--x; --y;
	assert(x >= 0 && x < t->n);
	assert(y >= 0 && y < t->m);
	va_list ap;
	memset(t->buf, 0, sizeof(t->buf));
	va_start(ap, fmt);
	vsnprintf(t->buf, sizeof(t->buf), fmt, ap);
	va_end(ap);
	t->cells[x][y] = strdup(t->buf);
}

void tab_display(tab_t *t) {
	int i, j, k;
	int *width = calloc(t->m, sizeof(int)); // cell width of each column.
	int table_w = 1;
	for(i = 0; i < t->m; ++i) {
		for(j = 0; j < t->n; ++j) {
			int len = t->cells[j][i] == NULL ?0 :(int)strlen(t->cells[j][i]);
			width[i] = width[i] > len ?width[i] :len; // space aligned and sufficient
		}
		table_w += width[i] + 2 + 1;
	}
	/* Output the table-header */
	if(t->header) {
		fprintf(stderr, "+");
		for(k = 0; k < table_w - 2; ++k) fprintf(stderr, "-");
		fprintf(stderr, "+\n");

		fprintf(stderr, "|");
		int spaces_n = table_w - 2 - (int)strlen(t->header);
		for(k = 0; k < spaces_n>>1; ++k) fprintf(stderr, " ");
		fprintf(stderr, "%s", t->header);
		for(k = 0; k < spaces_n - (spaces_n>>1); ++k) fprintf(stderr, " ");
		fprintf(stderr, "|\n");
	}
	/* Ouput table-cells  */
	for(i = 0; i < t->n; ++i) {
//		if(i == 0 || i == 1) {
			for(j = 0; j < t->m; ++j) {
				if(j == 0) fprintf(stderr, "+");
				for(k = 0; k < width[j] + 2; ++k) {
					fprintf(stderr, "-");
				}
				fprintf(stderr, "+");
			}
			fprintf(stderr, "\n");
//		}
		for(j = 0; j < t->m; ++j) {
			if(j == 0) fprintf(stderr, "|");
			int spaces_n = width[j] - (t->cells[i][j] == NULL ?0 :(int)strlen(t->cells[i][j]));
			for(k = 0; k < spaces_n>>1; ++k) fprintf(stderr, " ");
			if(t->cells[i][j]) fprintf(stderr, " %s ", t->cells[i][j]);
			else fprintf(stderr, "  ");
			for(k = 0; k < spaces_n - (spaces_n>>1); ++k) fprintf(stderr, " ");
			fprintf(stderr, "|");
		}
		fprintf(stderr, "\n");
		if(i == t->n - 1) {
			for(j = 0; j < t->m; ++j) {
				if(j == 0) fprintf(stderr, "+");
				for(k = 0; k < width[j] + 2; ++k) fprintf(stderr, "-");
				fprintf(stderr, "+");
			}
			fprintf(stderr, "\n");
		}
	}
	/* Output extra info */
	if(t->info) {
		fprintf(stderr, "%s", t->info);

		fprintf(stderr, "+");
		for(k = 0; k < table_w - 2; ++k) fprintf(stderr, "-");
		fprintf(stderr, "+\n");
	}
	free(width);
}

void tab_destroy(tab_t *t) {
	if(t->header) free(t->header);
	if(t->info) free(t->info);
	int i, j;
	for(i = 0; i < t->n; ++i) {
		for(j = 0; j < t->m; ++j) {
			free(t->cells[i][j]);
		}
		free(t->cells[i]);
	}
	free(t->cells);
}

void tab_test() {

	tab_t *t = tab_init(6, 7, "test-table", "11;\n22;\n");
	int i, j;
	for(i = 1; i <= t->n; ++i) {
		char str[20];
		for(j = 0; j < t->m; ++j) {
			str[j] = '*';
		}
		str[t->m] = '\0';
		for(j = 1; j <= t->m; ++j) {
			if (i == 1) {
				tab_fill_cell1(t, i, j, "%.3f", j * 123.0 + i);
			} else if (i == 2) {
				str[t->m - j + 1] = '\0';
				tab_fill_cell1(t, i, j, "%s", str);
			} else if (i == 3) {
				tab_fill_cell1(t, i, j, "%c", 'A' + j);
			} else if (i == 4) {
				tab_fill_cell1(t, i, j, "~%.6f X", j * 123.0 * i);
			} else if (i == 5) {
				tab_fill_cell1(t, i, j, "%d %%", j * 100 * i);
			} else if (i == 6) {
				tab_fill_cell1(t, i, j, "%d %%", j * 100 + i);
			}
		}
	}
	tab_display(t);
	tab_destroy(t);
	free(t);
}
