//
// Created by ixiaohu on 2019/10/29.
//

#ifndef TOY_1_1_PROGRESS_H
#define TOY_1_1_PROGRESS_H

#include <stdio.h>

typedef struct {
	char chr;		/*tip char*/
	char *title;	/*tip string*/
	int style;		/*progress style*/
	int max;		/*maximum value*/
	int last_pro;
	double offset;
	char *pro;
} progress_t;

#define PROGRESS_NUM_STYLE 0
#define PROGRESS_CHR_STYLE 1
#define PROGRESS_BGC_STYLE 2

void progress_init(progress_t *, char *, int, int);

void progress_show(progress_t *, double);

void progress_destroy(progress_t *);

#endif //TOY_1_1_PROGRESS_H
