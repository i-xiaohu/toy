//
// Created by ixiaohu on 2019/10/29.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "progress.h"

/**
 * initialize the progress bar.
 * @max = 0
 * @val = 0
 *
 * @param	style
 * @param	tip words.
 */
void progress_init(
		progress_t *bar, char *title, int max, int style)
{
	bar->chr = '#';
	bar->title = title;
	bar->style = style;
	bar->max = max;
	bar->last_pro = -1;
	bar->offset = 100 / (double)max;
	bar->pro = (char *) malloc((max + 1) * sizeof(char));
	if ( style == PROGRESS_BGC_STYLE )
		memset(bar->pro, 0x00, (max + 1) * sizeof(char));
	else {
		memset(bar->pro, 32, (max + 1) * sizeof(char));
		memset(bar->pro + max, 0x00, 1);
	}
}

void progress_show( progress_t *bar, double bit)
{
	int val = (int)(bit * bar->max);
	if(val == bar->last_pro) return ;
	bar->last_pro = val;
	switch ( bar->style )
	{
		case PROGRESS_NUM_STYLE:
			printf("\033[?25l\033[31m\033[1m%s%d%%\033[?25h\033[0m\r",
				   bar->title, (int)(bar->offset * val));
			fflush(stdout);
			break;
		case PROGRESS_CHR_STYLE:
			memset(bar->pro, '#', val * sizeof(char));
			printf("\033[?25l\033[31m\033[1m%s[%-s] %d%%\033[?25h\033[0m\r",
				   bar->title, bar->pro, (int)(bar->offset * val));
			fflush(stdout);
			break;
		case PROGRESS_BGC_STYLE:
			memset(bar->pro, 32, val * sizeof(char));
			printf("\033[?25l\033[31m\033[1m%s\033[41m %d%% %s\033[?25h\033[0m\r",
				   bar->title, (int)(bar->offset * val), bar->pro);
			fflush(stdout);
			break;
		default:
			break;
	}
}

//destroy the the progress bar.
void progress_destroy(progress_t *bar)
{
	free(bar->pro);
}