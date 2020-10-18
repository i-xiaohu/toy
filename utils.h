//
// Created by 63175 on 2019/7/31.
//

#ifndef TOY_UTILS_H
#define TOY_UTILS_H

#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>

static inline double cputime() {
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static inline double realtime() {
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

/* 高精度时间记录，time split专用 */
static inline double get_time() {
	struct timespec dida;
	clock_gettime(CLOCK_REALTIME, &dida);
	return (dida.tv_sec + dida.tv_nsec*1e-9);
}
#endif //TOY_UTILS_H
