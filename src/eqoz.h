#ifndef __RISM3D_EQOZ_H
#define __RISM3D_EQOZ_H

typedef void closure_t(int, const double *, const double *, double *, float *);
typedef void closure_c_t(int, const double *, const double *, double *);

#include <stdlib.h>

struct MEMORYINFO
{
	size_t c;
	size_t nr;
	size_t eq;
	size_t lsolv;
	size_t jx;
};
typedef struct MEMORYINFO meminfo_t;

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
#include "water.h"
#include "mol.h"

typedef struct {
	rism_t rism;
	solv_t solv;
	grid_t grid;
	closure_t *closure;
	/* private */
	double *solver_data;
	double *eq_data;
} eqoz_t;

typedef struct {
	int natv, lxvva;
	double *symc;
	float *dcdg;
	double *xvva; /* float *xvva; */
	unsigned int *indga;
	grid_t *grid;
	/* private */
	float *solver_data;
	float *jx_data;
} lneq_t;

void hnc(int n, const double *uuv, const double *tuv, double *cuv, float *f);
void hnc_c(int n, const double *uuv, const double *tuv, double *cuv);
void plhnc(int n, const double *uuv, const double *tuv, double *cuv, float *f);
void plhnc_c(int n, const double *uuv, const double *tuv, double *cuv);
int mk_eqoz(eqoz_t *eq, box_t *g, water_t *water, mol_t *mol,
            double ljcut, double ccut, int *spd, double th);
void rm_eqoz(eqoz_t *eq);
void eqoz(eqoz_t *eq, double *tuv, double *d, float *f);
int mk_lneq(lneq_t *ln, eqoz_t *eq);
void rm_lneq(lneq_t *ln);
void Jx(lneq_t *eq, float *x, float *r);
int nr(eqoz_t *eq, double *t, double *tol, int *maxit);

#endif //__RISM3D_EQOZ_H
