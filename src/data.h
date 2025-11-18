#ifndef __3D_RISM_DATA_H
#define __3D_RISM_DATA_H

#include <fftw3.h>

typedef struct
{
	fftw_plan fwd;
	fftw_plan bwd;
	int n, m, a, b, nk, nr, *nn;
	double *data;
} grid_d_t;

typedef struct
{
	fftwf_plan fwd;
	fftwf_plan bwd;
	int n, m, a, b, nk, nr, *nn;
	float *data;
} grid_f_t;

typedef struct {
    double *asymcr;
    double *asymck;
    double *asymhr;
    double *asymhk;
    double *huvk0;
    double *uuv;
} rism_t;

typedef struct {
    int natv;
    double *xvva;
    unsigned int *indga;
    int lxvva;
    double *charge;
    double *charge_sp;
} solv_t;

typedef void closure_t(int, double *, double *, double *, float *);

typedef struct {
    rism_t rism;
    solv_t solv;
    grid_d_t grid;
    grid_f_t lngr;
	closure_t *closure;
} eqoz_t;

typedef struct {
    int natv, lxvva;
    float *dcdg;
    double *xvva; /* float *xvva; */
    unsigned int *indga;
    grid_f_t *grid;
} lneq_t;

#endif //__3D_RISM_DATA_H
