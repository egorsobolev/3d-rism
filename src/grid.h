#ifndef __RISM3D_GRID_H
#define __RISM3D_GRID_H

#include <fft.h>

#define GRID_DIM 3

struct Box
{
	int n[3];
	double l[3];
	double s[3];
};
typedef struct Box box_t;

struct Grid
{
	size_t nr, nk, nk2, na;
	size_t n[3];
	size_t fft_shape[3];
	double l[3], s[3];
	double *v;
	double *v2;
	double *a;
	int *ia;
};
typedef struct Grid grid_t;


#include "mol.h"

void mkbox(const double *, const double *, const mol_t *, box_t *);
int ginit(const box_t *, grid_t *);
void rm_grid_wvec(grid_t *);

#endif //__RISM3D_GRID_H
