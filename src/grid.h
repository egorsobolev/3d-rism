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
	size_t nr, nk, nk2;
	int n[3];
	size_t fft_shape[3];
	double l[3], s[3], volume;
};
typedef struct Grid grid_t;

struct WaveVectors
{
	size_t nk2, na;
	double volume;
	double *v;
	double *v2;
	double *a;
	int *ia;
};
typedef struct WaveVectors wvec_t;


#include "mol.h"

void mkbox(const double *, const double *, const mol_t *, box_t *);
void grid_init(const box_t *, grid_t *);
int mk_wavevectors(grid_t *, wvec_t *);
void null_wavevectors(wvec_t *);
void rm_wavevectors(wvec_t *);

#endif //__RISM3D_GRID_H
