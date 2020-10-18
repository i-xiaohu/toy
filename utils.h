//
// Created by 63175 on 2019/7/28.
//

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <sys/resource.h>
#include <sys/time.h>
#include <string>
using namespace std;

inline double cputime() {
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

inline double realtime() {
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

inline int StringToInt(string a) {
	int res = 0;
	for(int i=0; i<a.length(); ++i){
		res *= 10;
		res += a[i]-'0';
	}
	return res;
}

#endif //SRC_UTILS_H
