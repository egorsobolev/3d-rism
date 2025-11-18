#ifndef __RISM3D_TDYN_H
#define __RISM3D_TDYN_H

#include "grid.h"
#include "water.h"

typedef double musc_t(const grid_t *, const water_t *, const double *, const double *, const double *);

double etot(const grid_t *, const water_t *, const double *, const double *, const double *);
double musc_hnc(const grid_t *, const water_t *, const double *, const double *, const double *);
double musc_plhnc(const grid_t *, const water_t *, const double *, const double *, const double *);
double mugf(const grid_t *, const water_t *, const double *, const double *, const double *);
double ichi(const grid_t *, const water_t *, const double *, const double *, const double *);

#endif //__RISM3D_TDYN_H
