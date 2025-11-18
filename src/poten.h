#ifndef __RISM3D_POTEN_H
#define __RISM3D_POTEN_H

#include "grid.h"
#include "water.h"
#include "mol.h"

void uljuv(const grid_t *g, const water_t *w, const mol_t *m, double rcut, double *u);
void ucolu(const grid_t *g, const water_t *w, const mol_t *m, const int *spd, double rcut, double *uuv);
void ucoluv(const water_t *w, int n, const double *u, double *uuv);

#endif //__RISM3D_POTEN_H
