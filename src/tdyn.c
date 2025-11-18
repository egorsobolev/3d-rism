#include "tdyn.h"

/* excess energy */
double etot(const grid_t *g, const water_t *w, const double *uuv, const double *huv, const double *cuv)
{
	double v, e, a;
	int k, i, j;
	v = g->l[0] * g->l[1] * g->l[2] / g->nr;
	e = 0.0;
	k = 0;
	for (j = 0; j < w->natom; ++j) {
		a = v * w->n[j] * w->m.rho;
		for (i = 0; i < g->nr; ++i) {
			e += a * (huv[k] + 1.0) * uuv[k];
			++k;
		}
	}
	return e;
}

/* excess chemical potential in Singer-Chandler theory for HNC closure */
double musc_hnc(const grid_t *g, const water_t *w, const double *uuv, const double *huv, const double *cuv)
{
	double v, e, a;
	int k, i, j;
	v = g->l[0] * g->l[1] * g->l[2] / g->nr;
	e = 0.0;
	k = 0;
	for (j = 0; j < w->natom; ++j) {
		a = v * w->n[j] * w->m.rho;
		for (i = 0; i < g->nr; ++i) {
			e += a * (0.5 * huv[k] * (huv[k] - cuv[k]) - cuv[k]);
			++k;
		}
	}
	return e;
}
/* excess chemical potential in Singer-Chandler theory for HNC closure */
double musc_plhnc(const grid_t *g, const water_t *w, const double *uuv, const double *huv, const double *cuv)
{
	double v, e, a;
	int k, i, j;
	v = g->l[0] * g->l[1] * g->l[2] / g->nr;
	e = 0.0;
	k = 0;
	for (j = 0; j < w->natom; ++j) {
		a = v * w->n[j] * w->m.rho;
		for (i = 0; i < g->nr; ++i) {
			e += a * (0.5 * huv[k] * ((huv[k] <= 0.0) * huv[k] - cuv[k]) - cuv[k]);
			++k;
		}
	}
	return e;
}
/* excess chemical potential in gauss fluctuation theory */
double mugf(const grid_t *g, const water_t *w, const double *uuv, const double *huv, const double *cuv)
{
	double v, e, a;
	int k, i, j;
	v = g->l[0] * g->l[1] * g->l[2] / g->nr;
	e = 0.0;
	k = 0;
	for (j = 0; j < w->natom; ++j) {
		a = -v * w->n[j] * w->m.rho;
		for (i = 0; i < g->nr; ++i) {
			e += a * (0.5 * huv[k] + 1.0) * cuv[k];
			++k;
		}
	}
	return e;
}

/* inverse compressibility */
double ichi(const grid_t *g, const water_t *w, const double *uuv, const double *huv, const double *cuv)
{
	double v, cuv0, a;
	int k, i, j;
	v = g->l[0] * g->l[1] * g->l[2] / g->nr;

	cuv0 = 0.0;
	k = 0;
	for (j = 0; j < w->natom; ++j) {
		a = v * w->n[j] * w->m.rho;
		for (i = 0; i < g->nr; ++i) {
			cuv0 += a * cuv[k];
			++k;
		}
	}

	return w->compc * (1.0 - cuv0) / w->m.rho;
}
