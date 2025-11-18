#ifndef __RISM3D_GRID_H
#define __RISM3D_GRID_H

#include <fftw3.h>

struct Box
{
	int n[3];
	double l[3];
	double s[3];
};
typedef struct Box box_t;

struct Grid
{
	int nr, nk, nk2, na;
	int n[3];
	double l[3], s[3];
	double *v;
	double *v2;
	double *a;
	int *ia;
};
typedef struct Grid grid_t;

struct GridJ
{
	fftw_plan fwd;
	fftw_plan bwd;
	int n, m, a, b, nk, nr, *nn;
	double *data;
};
typedef struct GridJ dgrid_t;

struct GridEQ
{
	fftwf_plan fwd;
	fftwf_plan bwd;
	int n, m, a, b, nk, nr, *nn;
	float *data;
};
typedef struct GridEQ fgrid_t;


#include "mol.h"

void mkbox(const double *, const double *, const mol_t *, box_t *);
int ginit(const box_t *, grid_t *);
int mkdgrid(int n, int *nn, dgrid_t *p);
void rmdgrid(dgrid_t *p);
int mkfgrid(int n, int *nn, fgrid_t *p);
void rmfgrid(fgrid_t *p);


#endif //__RISM3D_GRID_H
