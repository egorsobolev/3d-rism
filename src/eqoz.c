#include "eqoz.h"

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <cblas.h>

void hnc(int n, const double *uuv, const double *tuv, double *cuv, float *f)
{
	int i;
	double a;
	for (i = 0; i < n; ++i) {
		a = exp(tuv[i] - uuv[i]) - 1.0;
		f[i] = (float) a;
		cuv[i] = a - tuv[i];
	}
}
void hnc_c(int n, const double *uuv, const double *tuv, double *cuv)
{
	int i;
	for (i = 0; i < n; ++i)
		cuv[i] = exp(tuv[i] - uuv[i]) - 1.0 - tuv[i];
}

void plhnc(int n, const double *uuv, const double *tuv, double *cuv, float *f)
{
	int i;
	double a, b;
	for (i = 0; i < n; ++i) {
		a = tuv[i] - uuv[i];
		b = exp((a < 0.0) * a) - 1.0;
		cuv[i] = b + (a >= 0.0) * a - tuv[i];
		f[i] = (float) b;
	}
}
void plhnc_c(int n, const double *uuv, const double *tuv, double *cuv)
{
	int i;
	double a;
	for (i = 0; i < n; ++i) {
		a = tuv[i] - uuv[i];
		cuv[i] = exp((a < 0.0) * a) + (a >= 0.0) * a - 1.0 - tuv[i];
	}
}


void eqoz(dgrid_t *grid, rism_t *rism, solv_t *solv, closure_t *closure, double *tuv, double *d, float *f)
{
	int i, j, k, l, n, m;
	double *r, *t, *s, *src, *dst, *b, a;

	n = solv->natv * grid->nr;
	m = solv->natv * grid->nk;

	r = malloc(n * sizeof(double));
	t = malloc(m * sizeof(double));
	s = malloc(m * sizeof(double));

	(*closure)(n, rism->uuv, tuv, r, f);

	cblas_dcopy(n, r, 1, d, 1);
	/*
	 * d = r + rism.asymcr(:) * solv.charge;
	 */
	dst = d;
	for (j = 0; j < solv->natv; ++j) {
		cblas_daxpy(grid->nr, solv->charge[j], rism->asymcr, 1, dst, 1);
		dst += grid->nr;
	}

	/*
	 * t = fftw_r2c(rism.wisdom, r);
	 * t = t - rism.asymck(:) * solv.charge;
	 */
	src = d;
	dst = t;
	for (j = 0; j < solv->natv; j++) {
		b = grid->data;
		for (i = 0; i < grid->m; ++i) {
			/*memcpy(b, src, eq->p->b * sizeof(float));*/
			cblas_dscal(grid->b, 1.0/grid->nr, src, 1);
			cblas_dcopy(grid->b, src, 1, b, 1);
			b += grid->b;
			src += grid->b;
		}
		dgrid_fwd(grid);
		/*memcpy(dst, eq->p->data, eq->p->nk * sizeof(float));*/
		cblas_dcopy(grid->nk, grid->data, 1, dst, 1);
		cblas_daxpy(grid->nk, -solv->charge[j], rism->asymck, 1, dst, 1);
		cblas_dscal(grid->nk, solv->symc[j], dst, 1);
		dst += grid->nk;
	}

	/*
	 * iga=rism.indga(fix(((1:rism.ngk)+1)/2));
	 * s = reshape(sum(repmat(t,[1, 1, solv.natom]).*rism.xvva(iga,:,:),2),rism.ngk,solv.natom);
	 */
	memset(s, 0, m * sizeof(double));
	l = 0;
	for (j = 0; j < m; j += grid->nk) {
		/* side elemnts */
		for (k = 0; k < j; k += grid->nk) {
			for (i = 0; i < grid->nk; ++i) {
				a = solv->xvva[solv->indga[i / 2] + l];
				s[j + i] += a * t[k + i];
				s[k + i] += a * t[j + i];
			}
			l += solv->lxvva;
		}
		/* diagonal elements */
		for (i = 0; i < grid->nk; ++i)
			s[j + i] += t[j + i] * solv->xvva[solv->indga[i / 2] + l];
		l += solv->lxvva;
	}
	/*
	memset(s, 0, m * sizeof(double));
	dst = s;
	l = 0;
	for (j = 0; j < solv->natv; ++j) {
		src = t;
		for (k = 0; k < solv->natv; ++k) {
			for (i = 0; i < grid->nk; ++i)
				dst[i] += src[i] * solv->xvva[solv->indga[i / 2] + l];
			src += grid->nk;
			l += solv->lxvva;
		}
		dst += grid->nk;
	}
	*/

	/*
	 * s(1:2,:)=s(1:2,:)+rism.huvk0(1:2,:);
	 * s = s + rism.asymhk(:) * (solv.charge_sp / solv.dielconst);
	 */
	dst = s;
	src = rism->huvk0;
	for (j = 0; j < solv->natv; ++j) {
		cblas_dscal(grid->nk, 1.0 / solv->symc[j], dst, 1);
		dst[0] += src[0];
		dst[1] += src[1];

		cblas_daxpy(grid->nk, solv->charge_sp[j], rism->asymhk, 1, dst, 1);

		dst += grid->nk;
		src += 2;
	}

	/*
	 * d = fftw_c2r(rism.wisdom, s);
	 */
	dst = d;
	src = s;
	for (j = 0; j < solv->natv; j++) {
		/*memcpy(eq->p->data, src, eq->p->nk * sizeof(float));*/
		cblas_dcopy(grid->nk, src, 1, grid->data, 1);
		src += grid->nk;
		dgrid_bwd(grid);
		b = grid->data;
		for (i = 0; i < grid->m; ++i) {
			/*for (k = 0; k < eq->p->b; ++k)
				dst[k] += b[k];*/
			cblas_dcopy(grid->b, b, 1, dst, 1);
			b += grid->b;
			dst += grid->b;
		}
	}
	/*
	 * d = d - rism.asymhr * (solv.charge_sp / solv.dielconst);
	 * d = -r - tuv
	 */
	dst = d;
	src = r;
	for (j = 0; j < solv->natv; j++) {
		cblas_daxpy(grid->nr, -solv->charge_sp[j], rism->asymhr, 1, dst, 1);
		cblas_daxpy(grid->nr, -1.0, src, 1, dst, 1);
		cblas_daxpy(grid->nr, -1.0, tuv, 1, dst, 1);

		dst += grid->nr;
		src += grid->nr;
		tuv += grid->nr;
	}

	free(s);
	free(t);
	free(r);
}

/* + 2 * n * sizeof(float) */
void Jx(lneq_t *eq, float *x, float *r)
{
	float *d, *t, *src, *dst, *b, a;
	int i, j, k, l;
	int n, m;

	n = eq->natv * eq->grid->nr;
	m = eq->natv * eq->grid->nk;

	d = malloc(m * sizeof(float));
	t = malloc(m * sizeof(float));

	for (i = 0; i < n; i++) {
		r[i] = eq->dcdg[i] * x[i] / eq->grid->nr;
	}

	src = r;
	dst = d;
	for (j = 0; j < eq->natv; j++) {
		b = eq->grid->data;
		for (i = 0; i < eq->grid->m; ++i) {
			/*memcpy(b, src, eq->grid->b * sizeof(float));*/
			cblas_scopy(eq->grid->b,src,1,b,1);
			b += eq->grid->b;
			src += eq->grid->b;
		}

		fgrid_fwd(eq->grid);
		/*memcpy(dst, eq->grid->data, eq->grid->nk * sizeof(float));*/
		cblas_scopy(eq->grid->nk, eq->grid->data, 1, dst, 1);
		cblas_sscal(eq->grid->nk, (float) eq->symc[j], dst, 1);
		dst += eq->grid->nk;
	}

	memset(t, 0, m * sizeof(float));
	l = 0;
	for (j = 0; j < m; j += eq->grid->nk) {
		/* side elemnts */
		for (k = 0; k < j; k += eq->grid->nk) {
			for (i = 0; i < eq->grid->nk; ++i) {
				a = (float) eq->xvva[eq->indga[i / 2] + l];
				t[j + i] += a * d[k + i];
				t[k + i] += a * d[j + i];
			}
			l += eq->lxvva;
		}
		/* diagonal elements */
		for (i = 0; i < eq->grid->nk; ++i)
			t[j + i] += d[j + i] * (float) eq->xvva[eq->indga[i / 2] + l];
		l += eq->lxvva;
	}
	/*
	memset(t, 0, m * sizeof(float));
	dst = t;
	l = 0;
	for (j = 0; j < eq->natv; ++j) {
		src = d;
		for (k = 0; k < eq->natv; ++k) {
			for (i = 0; i < eq->grid->nk; ++i)
				dst[i] += src[i] * (float) eq->xvva[eq->indga[i / 2] + l];
			src += eq->grid->nk;
			l += eq->lxvva;
		}
		dst += eq->grid->nk;
	}
	*/
	cblas_sscal(n, -1.0f * eq->grid->nr, r, 1);
	cblas_saxpy(n, -1.0f, x, 1, r, 1);
	dst = r;
	src = t;
	for (j = 0; j < eq->natv; j++) {
		cblas_sscal(eq->grid->nk, 1.0f / (float) eq->symc[j], src, 1);
		/*memcpy(eq->grid->data, src, eq->grid->nk * sizeof(float));*/
		cblas_scopy(eq->grid->nk, src, 1, eq->grid->data, 1);
		src += eq->grid->nk;
		fgrid_bwd(eq->grid);
		b = eq->grid->data;
		for (i = 0; i < eq->grid->m; ++i) {
			/*for (k = 0; k < eq->grid->b; ++k)
				dst[k] += b[k];*/
			cblas_saxpy(eq->grid->b,1.0f,b,1,dst,1);
			b += eq->grid->b;
			dst += eq->grid->b;
		}
	}

	free(t);
	free(d);
}
