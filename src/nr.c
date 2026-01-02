#include <memory.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "eqoz.h"
#include "utils.h"

#define print  printf
#define flush  fflush(stdin)

#define MAX_ARM  (5)
#define SCALE_ARM  (0.618033988749895)


/* + 5 * n * sizeof(float)*/
int bicgstab(lneq_t *eq, int N, const float *b, float *x, float *tol, int *it)
{
	float *r, *h, *p, *v, *s;
	float rTh, rTr, sTr, hTv;
	float sumb2, sumr2, sums2;
	int its, maxit, more, exitcode;

	/*
	 r, h, p, v, s: float[N]
	 */
	r = eq->solver_data;
	h = r + N;
	p = h + N;
	v = p + N;
	s = v + N;

	sumb2 = 0.0f;
	rTh = 0.0f;
	sumr2 = 0.0f;
	#pragma omp parallel
	{
		int i, maxit, its, more;
		float normb, normr, norms;
		float err, omega, beta, alpha;

		maxit = *it * 2;
		its = 0;
		more = 1;
		
		#pragma omp for reduction(+:sumb2)
		for (i = 0; i < N; i++)
			sumb2 += b[i] * b[i];
		normb = sqrtf(sumb2);
		if (normb == 0.0f)
			normb = 1.0f;
		err = *tol * normb;
		Jx(eq, x, r);
		++its;

		#pragma omp for reduction(+:rTh, sumr2)
		for (i = 0; i < N; i++) {
			r[i] -= b[i];
			p[i] = h[i] = r[i];
			rTh += r[i] * h[i];
			sumr2 += r[i] * r[i];
		}
		normr = sqrtf(sumr2);
		if (normr <= err) {
			#pragma omp single
			{
				*tol = normr / normb;
				*it = its;
				exitcode = 0;
			}
			more = 0;
		}
		while (more && its < maxit) {
			Jx(eq,p,v);  /* Jx(eq, p, v) */
			++its;
			#pragma omp single nowait
			hTv = 0.0f;
			#pragma omp for reduction(+:hTv)
			for (i = 0; i < N; i++)
				hTv += h[i] * v[i];
			alpha = rTh / hTv;
			#pragma omp single nowait
			sums2 = 0.0f;
			#pragma omp for reduction(+:sums2)
			for (i = 0; i < N; i++) {
				s[i] = r[i] - alpha * v[i];
				sums2 += s[i];
			}
			norms = sqrtf(sums2);
			if (norms <= err) {
				#pragma omp for
				for (i = 0; i < N; i++)
					x[i] -= alpha * p[i];
				#pragma omp single
				{
					*tol = norms / normb;
					*it = its;
					exitcode = 0;
				}
				more = 0;
				break;
			}
			Jx(eq,s,r);
			++its;
			#pragma omp single nowait
			sTr = 0.0f;
			#pragma omp single nowait
			rTr = 0.0f;
			#pragma omp for reduction(+:sTr, rTr)
			for (i = 0; i < N; i++) {
				sTr += s[i] * r[i];
				rTr += r[i] * r[i];
			}
			if ( fabs(sTr)<1e-40 || fabs(rTr)<1e-40 )
				omega = 0.;
			else
				omega = sTr/rTr;
			#pragma omp single nowait
			sumr2 = 0.0f;
			#pragma omp for reduction(+:sumr2)
			for (i = 0; i < N; i++) {
				x[i] -= alpha * p[i] + omega * s[i];
				r[i] = s[i] -omega * r[i];
				sumr2 += r[i] * r[i];
			}
			normr = sqrtf(sumr2);
			if (normr <= err) {
				#pragma omp single
				{
					*tol = normr / normb;
					*it = its;
					exitcode = 0;
				}
				more = 0;
				break;
			}
			if (omega == 0.0f) {
				#pragma omp single
				{
					*tol = normr / normb;
					*it = its;
					exitcode = 2;
				}
				more = 0;
				break;
			}
			beta=(alpha/omega)/rTh;
			#pragma omp barrier
			#pragma omp single nowait
			rTh = 0.0f;
			#pragma omp for reduction(+:rTh)
			for (i = 0; i < N; i++)
				rTh += r[i] * h[i];
			if (rTh == 0.0f) {
				#pragma omp single
				{
					*tol = normr / normb;
					*it = its;
					exitcode = 3;
				}
				more = 0;
				break;
			}
			beta*=rTh;
			#pragma omp for
			for (i = 0; i < N; i++)
				p[i] = beta * p[i] + r[i] - beta * omega * v[i];
		}
		if (more) {
			#pragma omp single
			{
				*tol = rTr / normb;
				*it = its;
				exitcode = 1;
			}
		}
	}
	return exitcode; /* no convergence */
}

/* Eisenstat, Walker */
typedef struct
{
	double etaM;
	double eta;
	double gamma;
	double prec;
	double normdT;
} fterm_ew_t;

void ftew_setdef(fterm_ew_t *f, double tol)
{
	f->etaM = 0.9;
	f->eta = sqrt(f->etaM / 0.9);
	f->gamma = 0.9;
	f->prec = tol;
	f->normdT = 1.0;
}

double ftew_next(fterm_ew_t *f, double normd)
{
	double etaT;
	etaT = f->gamma * f->eta * f->eta;
	f->eta = normd / f->normdT;
	f->eta = f->gamma * f->eta * f->eta;
	if (etaT > .1 && etaT > f->eta)
		f->eta = etaT;
	if (f->etaM < f->eta)
		f->eta = f->etaM;
	etaT = f->prec / normd;
	if (etaT > f->eta)
		f->eta = etaT;
	f->normdT = normd;

	return f->eta;
}

/* 9 * n * sizeof(float) */
int nr(eqoz_t *eq, double *t, double *tol, int *maxit)
{
	double *d, *s, *tN, norms, eta, sumd2, sqrt_n;
	float *b, *x, eps;
	int n, m, i, j, flag, nlit, maxlit, fail, narm;
	double lntm, nrtm, lntmi, nrtm_cpu, err;
	fterm_ew_t fterm;

	nrtm_cpu = clock() / (double) CLOCKS_PER_SEC;
	nrtm = walltime() * 1e-6;
	lntm = .0;

	ftew_setdef(&fterm, *tol);
	maxlit = 1000;

	n = eq->solv.natv * eq->grid.nr;
	sqrt_n = sqrt(n);

	/*
	 b, x: float[n]
	 d, s, tN: double[n]
	 */
	d = eq->solver_data;
	s = d + n;
	tN = s + n;
	b = (float *) (tN + n);
	x = b + n;
	
	fail = 0;
	sumd2 = 0.0;
	#pragma omp parallel
	{
		int k;

		eqoz(eq, t, d, eq->lneq.dcdg);

		#pragma omp for reduction(+:sumd2)
		for (k = 0; k < n; k++)
			sumd2 += d[k] * d[k];
		#pragma omp single nowait
		err = sqrt(sumd2) / sqrt_n;
		
		#pragma omp for
		for (k = 0; k < n; ++k) {
			b[k] = (float) d[k];
			x[k] = .0f;
		}
	}

	print("||Z(t0)|| = %.1e\n", err);
	print("    #     Nit F      t,s M   eta   ||Jx-Z||       ||Z||\n");
	flush;
	i = 0;
	while (i <= *maxit && err > *tol) {

		nlit = maxlit;
		eta = ftew_next(&fterm, err);
		lntmi = walltime() * 1e-6;
		eps = (float) eta;
		flag = bicgstab(&eq->lneq, n, b, x, &eps, &nlit);
		lntmi = walltime() * 1e-6 - lntmi;
		lntm += lntmi;

		if (flag) {
			/* no convergence of linear solver */
			return 1;
		}
		norms = (double) eps;
		
		#pragma omp parallel
		{
			int k, j;
			double errN, lambda;

			#pragma omp for
			for (k = 0; k < n; k++) {
				s[k] = (double) -x[k];
				tN[k] = t[k] + s[k];
			}

			eqoz(eq, tN, d, eq->lneq.dcdg);
			#pragma single nowait
			sumd2 = 0.0;
			#pragma omp for reduction(+:sumd2)
			for (k = 0; k < n; k++)
				sumd2 += d[k] * d[k];
			errN = sqrt(sumd2) / sqrt_n;

			j = 0;
			lambda = 1.0;
			while (j < MAX_ARM && errN > err) {
				lambda *= SCALE_ARM;
				#pragma omp for
				for (k = 0; k < n; k++)
					tN[k] = t[k] + lambda * s[k];

				eqoz(eq, tN, d, eq->lneq.dcdg);

				#pragma omp single nowait
				sumd2 = 0.0;
				#pragma omp for reduction(+:sumd2)
				for (k = 0; k < n; k++)
					sumd2 += d[k] * d[k];
				errN = sqrt(sumd2) / sqrt_n;
				++j;
			}

			if (j >= MAX_ARM) {
				/* no convergence of armigo*/
				#pragma omp single nowait
				fail = 1;
			} else {
				#pragma omp for
				for (k = 0; k < n; k++) {
					t[k] = tN[k];
					b[k] = (float) d[k];
					x[k] = .0f;
				}
				#pragma omp single nowait
				err = errN;
				#pragma omp single nowait
				narm = j;
			}
		}
		if (fail) {
			/* no convergence of armigo*/
			return 2;
		}
		print("%5d %7d %1d %8.1lf %1d %5.2lf %10.6lf %11.6lf\n",
		      i+1, nlit, flag, lntmi, narm, eta, norms, err);
		flush;
		++i;
	}
	*maxit = i;
	*tol = err;

	nrtm_cpu = clock() / (double) CLOCKS_PER_SEC - nrtm_cpu;
	nrtm = walltime() * 1e-6 - nrtm;
	print("NR: %.1lfs, BiCGStab: %.1lfs, CPU usage: %.1lf\n", nrtm, lntm, nrtm_cpu / nrtm);
	print("||Z(tn)|| = %.1e\n", *tol);

	return 0; /* OK */
}
