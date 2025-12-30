#ifndef __3D_RISM_DATA_H
#define __3D_RISM_DATA_H

struct Grid
{
	int nr, nk, nk2, na;
	int n[3];
	size_t fft_shape[3];
	double l[3], s[3];
	double *v;
	double *v2;
	double *a;
	int *ia;
};
typedef struct Grid grid_t;

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
	grid_t grid;
	closure_t *closure;
} eqoz_t;

typedef struct {
	int natv, lxvva;
	float *dcdg;
	double *xvva; /* float *xvva; */
	unsigned int *indga;
	grid_t *grid;
} lneq_t;

#endif //__3D_RISM_DATA_H
