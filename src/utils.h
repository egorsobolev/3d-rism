#ifndef __RISM3D_UTILS_H
#define __RISM3D_UTILS_H

#define PARABOLA	0
#define FDERIV		1
#define SDERIV		2

int cubicuni(int, double *, double *, int, double, int, double);
void resample(int, const double *, int, const double *, double *);
void quicksort(int, double *, int *);

#endif //__RISM3D_UTILS_H
