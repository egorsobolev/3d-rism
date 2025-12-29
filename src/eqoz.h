#ifndef __RISM3D_EQOZ_H
#define __RISM3D_EQOZ_H

typedef void closure_t(int, const double *, const double *, double *, float *);
typedef void closure_c_t(int, const double *, const double *, double *);

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
	double *symc;
} solv_t;

#include "grid.h"

typedef struct {
	rism_t rism;
	solv_t solv;
	dgrid_t grid;
	fgrid_t lngr;
	closure_t *closure;
} eqoz_t;

typedef struct {
	int natv, lxvva;
	double *symc;
	float *dcdg;
	double *xvva; /* float *xvva; */
	unsigned int *indga;
	fgrid_t *grid;
} lneq_t;

void hnc(int n, const double *uuv, const double *tuv, double *cuv, float *f);
void hnc_c(int n, const double *uuv, const double *tuv, double *cuv);
void plhnc(int n, const double *uuv, const double *tuv, double *cuv, float *f);
void plhnc_c(int n, const double *uuv, const double *tuv, double *cuv);
void eqoz(dgrid_t *grid, rism_t *rism, solv_t *solv, closure_t *closure, double *tuv, double *d, float *f);
void Jx(lneq_t *eq, float *x, float *r);
int nr(eqoz_t *eq, double *t, double *tol, int *maxit);

#endif //__RISM3D_EQOZ_H
