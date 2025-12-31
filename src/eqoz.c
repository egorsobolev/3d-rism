#include "eqoz.h"

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <cblas.h>

#include "utils.h"
#include "asymp.h"
#include "poten.h"


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


void null_solv(solv_t *solv)
{
	solv->charge = NULL;
}

int mk_solv(solv_t *solv, wvec_t *wvec, water_t *water)
{
	int err;
	walltime_t t0;

	solv->natv = water->natom;

	/*
	 charge: double[natv]
	 charge_sp: double[natv]
	 symc: double[natv]
	 xvva: double[3 * na]
	 indga: unsigned int[nk]
	 */
	solv->charge = (double *) malloc(3 * (solv->natv + wvec->na) * sizeof(double) +
	                                 wvec->nk2 * sizeof(int));
	if (!solv->charge)
		return -1;
	solv->charge_sp = solv->charge + solv->natv;
	solv->symc = solv->charge_sp + solv->natv;
	solv->xvva = solv->symc + solv->natv;
	solv->indga = (int *) (solv->xvva + 3 * wvec->na);

	/* NOTE: only water */
	solv->charge[0] = water->m.q_h;
	solv->charge[1] = -2.0 * water->m.q_h;
	solv->charge_sp[0] = 0.0;
	solv->charge_sp[1] = 0.0;
	solv->symc[0] = sqrt(water->n[0]);
	solv->symc[1] = sqrt(water->n[1]);

	solv->lxvva = wvec->na;
	memcpy(solv->indga, wvec->ia, wvec->nk2 * sizeof(int));

	t0 = walltime();
	err = mkxvva(wvec, water, solv->xvva);
	if (err == 1) {
		printf(" MKXVVA: solvent functions are too short\n");
		return -2;
	} else if (err == 2) {
		printf(" MKXVVA: insufficient memory\n");
		return -3;
	}
	printf(" MKXVVA: Elapse %.2lf second(s)\n", (walltime() - t0) * 1e-6);

	return 0;
}


void rm_solv(solv_t *solv)
{
	free(solv->charge);
}


void null_rism(rism_t * rism)
{
	rism->uuv = NULL;
}


int mk_rism(rism_t *rism, grid_t *g, wvec_t *wvec, water_t *water, mol_t *mol,
            double ljcut, double ccut, int *spd, double th)
{
	size_t size, n, m;
	double *u0;
	walltime_t t0;

	n = water->natom * g->nr;
	//m = water->natom * g->nk;
	/*
	 uuv: double[natv * nr]
	 asymcr: double[nr]
	 asymhr: double[nr]
	 asymck: double[nk]
	 asymhk: double[nk]
	 huvk0: double[2 * natv]
	 */

	u0 = (double *) malloc(g->nr * sizeof(double));
	if (!u0) return -1;

	size = n + 2 * (g->nr + g->nk + water->natom);
	rism->uuv = (double *) malloc(size * sizeof(double));
	if (!rism->uuv) {
		free(u0);
		return -1;
	}

	rism->asymcr = rism->uuv + n;
	rism->asymhr = rism->asymcr + g->nr;
	rism->asymck = rism->asymhr + g->nr;
	rism->asymhk = rism->asymck + g->nk;
	rism->huvk0 = rism->asymhk + g->nk;

	t0 = walltime();
	memset(rism->uuv, 0, n * sizeof(double));
	uljuv(g, water, mol, ljcut, rism->uuv);
	printf(" ULJUV : Elapse %.2lf second(s)\n", (walltime() - t0) * 1e-6);

	t0 = walltime();
	ucolu(g, water, mol, spd, ccut, u0);
	ucoluv(water, g->nr, u0, rism->uuv);
	printf(" UCOLUV: Elapse %.2lf second(s)\n", (walltime() - t0) * 1e-6);

	t0 = walltime();
	asympr(g, water, mol, u0, th, rism->asymcr, rism->asymhr);
	printf(" ASYMPR: Elapse %.2lf second(s)\n", (walltime() - t0) * 1e-6);

	t0 = walltime();
	asympk(wvec, water, mol, th, rism->asymck, rism->asymhk, rism->huvk0);
	printf(" ASYMPK: Elapse %.2lf second(s)\n", (walltime() - t0) * 1e-6);
	printf("\n");

	free(u0);
	return 0;
}


void rm_rism(rism_t *rism)
{
	free(rism->uuv);
}


void null_lneq(lneq_t *ln)
{
	ln->dcdg = NULL;
}

int mk_lneq(lneq_t *ln, eqoz_t *eq)
{
	size_t natv, nr, nk, n, m;

	natv = eq->solv.natv;
	nr = eq->grid.nr;
	nk = eq->grid.nk;
	n = nr * natv;
	m = nk * natv;

	ln->natv = natv;
	ln->symc = eq->solv.symc;
	ln->lxvva = eq->solv.lxvva;
	ln->xvva = eq->solv.xvva;
	ln->indga = eq->solv.indga;
	ln->grid = &eq->grid;

	/*
	 dcdg: float[n]
	 jx_data: float[2 * m]
	 solver_data: float[5 * n]
	 */
	ln->dcdg = (float *) malloc((6*n + 2*m) * sizeof(float));
	if (!ln->dcdg)
		return -1;

	ln->jx_data = ln->dcdg + n;
	ln->solver_data = ln->jx_data + 2*m;

	return 0;
}


void rm_lneq(lneq_t *ln)
{
	free(ln->dcdg);
}


void print_mem_usage(grid_t *g, wvec_t *wvec, water_t *w, mol_t *m)
{
	meminfo_t ms, mn;

	ms.c = (3 * wvec->na + 3 * g->nk + 2 * g->nr + w->natom * g->nr + 6 * m->natom) * sizeof(double) + 
	       g->nk * (sizeof(int) + sizeof(float));
	mn.c = ms.c + (wvec->na + 4 * g->nk + 3 * w->ngrid) * sizeof(double);

	ms.nr = (4 * sizeof(double) + 3 * sizeof(float)) * g->nr * w->natom;
	mn.nr = ms.nr + ms.c;

	ms.eq = (2 * g->nk + g->nr) * sizeof(double) * w->natom;
	mn.eq = mn.nr + ms.eq;

	ms.lsolv = 5 * g->nr * sizeof(float) * w->natom;
	mn.lsolv = mn.nr + ms.lsolv;

	ms.jx = 2 * g->nk * sizeof(float) * w->natom;
	mn.jx = mn.lsolv + ms.jx;

	printf(" MEMORY            self, B       net, B\n");
	printf("  constants   %12ld %12ld\n", ms.c, mn.c);
	printf("   NR         %12ld %12ld\n", ms.nr, mn.nr);
	printf("    EQ        %12ld %12ld\n", ms.eq, mn.eq);
	printf("    BiCGStab  %12ld %12ld\n", ms.lsolv, mn.lsolv);
	printf("     Jx       %12ld %12ld\n", ms.jx, mn.jx);
	printf("\n");
}


int mk_eqoz(eqoz_t *eq, box_t *box, water_t *water, mol_t *mol,
            double ljcut, double ccut, int *spd, double th)
{
	size_t natv, nk, nr, m, n;
	size_t solver_data_size;
	size_t eq_data_size;
	wvec_t wvec;
	walltime_t t0;

	null_wavevectors(&wvec);
	null_solv(&eq->solv);
	null_rism(&eq->rism);
	null_lneq(&eq->lneq);
	eq->eq_data = NULL;

	grid_init(box, &eq->grid);

	t0 = walltime();
	if (mk_wavevectors(&eq->grid, &wvec)) {
		printf(" MK_WAVEVECTORS : insufficient memory\n");
		return -2;
	}
	printf(" && %ld, %lf\n", wvec.na, wvec.v2[wvec.na - 2]);
	t0 = walltime() - t0;
	print_mem_usage(&eq->grid, &wvec, water, mol);
	printf(" GINIT : Elapse %.2lf second(s)\n", t0 * 1e-6);

	/* TODO: catch exceptions */
	if (mk_solv(&eq->solv, &wvec, water)) {
		printf("MK_SOLV: insufficient memory\n");
		rm_wavevectors(&wvec);
		return -3;
	}
	if (mk_rism(&eq->rism, &eq->grid, &wvec, water, mol, ljcut, ccut, spd, th)) {
		printf("MK_RISM: insufficient memory\n");
		rm_wavevectors(&wvec);
		return -4;
	}

	natv = eq->solv.natv;
	nr = eq->grid.nr;
	nk = eq->grid.nk;
	n = nr * natv;
	m = nk * natv;

	/*
	 solver_data: float[2 * n] & double[3 * n]
	 eq_data: double[n + 2 * m]
	 */
	solver_data_size = n * (2 * sizeof(float) + 3 * sizeof(double));
	eq_data_size = (n + 2 * m) * sizeof(double);
	eq->eq_data = (double *) malloc(solver_data_size + eq_data_size);
	if (!eq->eq_data) {
		printf("MK_EQOZ: insufficient memory\n");
		rm_wavevectors(&wvec);
		return -1;
	}
	eq->solver_data = eq->eq_data + eq_data_size / sizeof(double);

	if (mk_lneq(&eq->lneq, eq)) {
		printf("MK_LNEQ: insufficient memory\n");
		rm_wavevectors(&wvec);
		return -1;
	}

	rm_wavevectors(&wvec);
	return 0;
}


void rm_eqoz(eqoz_t *eq)
{
	rm_lneq(&eq->lneq);
	free(eq->eq_data);
	rm_rism(&eq->rism);
	rm_solv(&eq->solv);
}


void eqoz(eqoz_t *eq, double *tuv, double *d, float *f)
{
	int i, j, k, l, n, m;
	double *r, *t, *s, *src, *dst, *b, a;
	size_t nk, nr, natv;

	nk = eq->grid.nk;
	nr = eq->grid.nr;
	natv = eq->solv.natv;

	n = natv * nr;
	m = natv * nk;

	r = eq->eq_data;
	t = r + n;
	s = t + m;

	(*eq->closure)(n, eq->rism.uuv, tuv, r, f);

	cblas_dcopy(n, r, 1, d, 1);
	/*
	 * d = r + rism.asymcr(:) * solv.charge;
	 */
	dst = d;
	for (j = 0; j < natv; ++j) {
		cblas_daxpy(nr, eq->solv.charge[j], eq->rism.asymcr, 1, dst, 1);
		dst += nr;
	}

	/*
	 * t = fftw_r2c(rism.wisdom, r);
	 * t = t - rism.asymck(:) * solv.charge;
	 */
	src = d;
	dst = t;
	for (j = 0; j < natv; j++) {
		cblas_dscal(nr, 1.0/nr, src, 1);
		fft_r2c(3, eq->grid.fft_shape, src, (complex_double *) dst, 1.0);
		cblas_daxpy(nk, -eq->solv.charge[j], eq->rism.asymck, 1, dst, 1);
		cblas_dscal(nk, eq->solv.symc[j], dst, 1);

		src += nr;
		dst += nk;
	}

	/*
	 * iga=rism.indga(fix(((1:rism.ngk)+1)/2));
	 * s = reshape(sum(repmat(t,[1, 1, solv.natom]).*rism.xvva(iga,:,:),2),rism.ngk,solv.natom);
	 */
	memset(s, 0, m * sizeof(double));
	l = 0;
	for (j = 0; j < m; j += nk) {
		/* side elemnts */
		for (k = 0; k < j; k += nk) {
			for (i = 0; i < nk; ++i) {
				a = eq->solv.xvva[eq->solv.indga[i / 2] + l];
				s[j + i] += a * t[k + i];
				s[k + i] += a * t[j + i];
			}
			l += eq->solv.lxvva;
		}
		/* diagonal elements */
		for (i = 0; i < nk; ++i)
			s[j + i] += t[j + i] * eq->solv.xvva[eq->solv.indga[i / 2] + l];
		l += eq->solv.lxvva;
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
	src = eq->rism.huvk0;
	for (j = 0; j < natv; ++j) {
		cblas_dscal(nk, 1.0 / eq->solv.symc[j], dst, 1);
		dst[0] += src[0];
		dst[1] += src[1];

		cblas_daxpy(nk, eq->solv.charge_sp[j], eq->rism.asymhk, 1, dst, 1);

		dst += nk;
		src += 2;
	}

	/*
	 * d = fftw_c2r(rism.wisdom, s);
	 */
	dst = d;
	src = s;
	for (j = 0; j < natv; j++) {
		fft_c2r(3, eq->grid.fft_shape, (complex_double *) src, dst, 1.0);
		src += nk;
		dst += nr;
	}
	/*
	 * d = d - rism.asymhr * (solv.charge_sp / solv.dielconst);
	 * d = -r - tuv
	 */
	dst = d;
	src = r;
	for (j = 0; j < natv; j++) {
		cblas_daxpy(nr, -eq->solv.charge_sp[j], eq->rism.asymhr, 1, dst, 1);
		cblas_daxpy(nr, -1.0, src, 1, dst, 1);
		cblas_daxpy(nr, -1.0, tuv, 1, dst, 1);

		dst += nr;
		src += nr;
		tuv += nr;
	}
}


/* + 2 * n * sizeof(float) */
void Jx(lneq_t *eq, float *x, float *r)
{
	float *d, *t, *src, *dst, *b, a;
	int i, j, k, l;
	int n, m;
	size_t nk, nr;

	nk = eq->grid->nk;
	nr = eq->grid->nr;
	n = eq->natv * nr;
	m = eq->natv * nk;

	d = eq->jx_data;
	t = d + m;

	for (i = 0; i < n; i++) {
		r[i] = eq->dcdg[i] * x[i] / nr;
	}

	src = r;
	dst = d;
	for (j = 0; j < eq->natv; j++) {
		fftf_r2c(3, eq->grid->fft_shape, src, (complex_float *) dst, 1.0f);
		cblas_sscal(nk, (float) eq->symc[j], dst, 1);

		src += nr;
		dst += nk;
	}

	memset(t, 0, m * sizeof(float));
	l = 0;
	for (j = 0; j < m; j += nk) {
		/* side elemnts */
		for (k = 0; k < j; k += nk) {
			for (i = 0; i < nk; ++i) {
				a = (float) eq->xvva[eq->indga[i / 2] + l];
				t[j + i] += a * d[k + i];
				t[k + i] += a * d[j + i];
			}
			l += eq->lxvva;
		}
		/* diagonal elements */
		for (i = 0; i < nk; ++i)
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
	cblas_sscal(n, -1.0f * nr, r, 1);
	cblas_saxpy(n, -1.0f, x, 1, r, 1);
	dst = r;
	src = t;
	for (j = 0; j < eq->natv; j++) {
		cblas_sscal(nk, 1.0f / (float) eq->symc[j], src, 1);
		fftf_c2r(3, eq->grid->fft_shape, (complex_float *) src, d, 1.0f);
		cblas_saxpy(nr,1.0f,d,1,dst,1);

		src += nk;
		dst += nr;
	}
}
