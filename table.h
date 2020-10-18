//
// Created by 63175 on 2019/12/1.
//

#ifndef TOY_1_1_TABLE_H
#define TOY_1_1_TABLE_H

typedef struct {
	int n, m; // n rows, m columns
	char *header, *info;
	char ***cells;
	char buf[256]; // for sprintf
} tab_t;

#ifdef __cplusplus
extern "C" {
#endif
	tab_t *tab_init(int n, int m, char *header, char *info);
	void tab_test();
	void tab_fill_cell1(tab_t *t, int x, int y, const char *fmt, ...);
	void tab_display(tab_t *t);
	void tab_destroy(tab_t *t);

#ifdef __cplusplus
}
#endif

#endif //TOY_1_1_TABLE_H
