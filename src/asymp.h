#ifndef __RISM3D_ASYMP_H
#define __RISM3D_ASYMP_H

#include "grid.h"
#include "water.h"
#include "mol.h"

void asympr(const grid_t *g, const water_t *w, const mol_t *m, const double *u, double th, double *asympcr, double *asymphr);
void asympk(const grid_t *g, const water_t *w, const mol_t *m, double th, double *asympck, double *asymphk, double *huvk0);

#endif //__RISM3D_ASYMP_H
