#ifndef __RISM3D_UTILS_H
#define __RISM3D_UTILS_H

#include <stdint.h>

#define PARABOLA	0
#define FDERIV		1
#define SDERIV		2

int cubicuni(int, double *, double *, int, double, int, double);
void resample(int, const double *, int, const double *, double *);
void quicksort(int, double *, int *);

typedef int64_t walltime_t;
walltime_t walltime();

#endif //__RISM3D_UTILS_H
