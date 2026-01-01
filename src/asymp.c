#include "asymp.h"

#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cblas.h>

#define RCOR		(0.002)
#define M_4PI		(4.0 * M_PI)
#define M_1_2PI2	(1.0 / (2.0 * M_PI * M_PI))
#define MAX(a,b)	(a > b ? a : b);

#ifndef erfc
double erfc(double x)
{
	double t, z, ans;
	z = fabs(x);
	t = 1.0 / (1.0 + 0.5 * z);

	ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	     t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	     t*(-0.82215223+t*0.17087277)))))))));

	return x >= 0.0 ? ans : 2.0-ans;
}
#endif

void vec_sincos(size_t n, const double *x, double *s, double *c)
{
	for (size_t i = 0; i < n; ++i) {
		s[i] = sin(x[i]);
		c[i] = cos(x[i]);
	}
}

void asympr(const grid_t *g, const water_t *w, const mol_t *m, const double *u, double th, double *asympcr, double *asymphr)
{
	int i, l, k;
	int x, x0, xn, y, y0, yn, z, z0, zn;
	double smear, xappa, a, b, d, q_r;
	double rcut, rcut2, rs, ccuthr, ccutcr, thhr, thcr, casympcr, casymphr;
	double r2, ra, dz2, dz2dy2;
	smear = 1.0;
	xappa = 0.0;

	a = 0.5 * xappa * smear;
	b = exp(-a * a);
	thhr = b - th;
	thcr = 1.0 - th;
	ccuthr = 0.0;
	ccutcr = 0.0;
	rcut = 0.0;
	do {
		rcut += 1.0;

		d = exp(-xappa * rcut);
		rs = rcut / smear;

		ccuthr = 0.5 * b * (d * erfc(a - rs) - erfc(a + rs) / d);
		ccutcr = (1.0 - erfc(rs));
	} while (thhr > ccuthr && thcr > ccutcr);
	casympcr = M_2_SQRTPI / smear;
	casymphr = b * (casympcr * b - xappa * erfc(a));

	rcut2 = rcut * rcut;
	
	for (k = 0; k < g->nr; ++k) {
		asympcr[k] = ccutcr * u[k];
		asymphr[k] = ccuthr * u[k];
	}

	for (i = 0; i < m->natom; ++i) {

		z0 = (int) ceil((m->z[i] - rcut) / g->s[2]);
		if (z0 < 0) z0 = 0;
		zn = (int) floor((m->z[i] + rcut) / g->s[2]);
		if (zn >= g->n[2]) zn = g->n[2] - 1;

		for (z = z0; z <= zn; ++z) {
			d = z * g->s[2] - m->z[i];
			dz2 = d * d;

			d = sqrt(rcut2 - dz2);
			y0 = (int) ceil((m->y[i] - d) / g->s[1]);
			if (y0 < 0) y0 = 0;
			yn = (int) floor((m->y[i] + d) / g->s[1]);
			if (yn >= g->n[1]) yn = g->n[1] - 1;

			l = z * g->n[1];
			for (y = y0; y <= yn; ++y) {
				d = y * g->s[1] - m->y[i];
				dz2dy2 = dz2 + d * d;

				d = sqrt(rcut2 - dz2dy2);
				x0 = (int) ceil((m->x[i] - d) / g->s[0]);
				if (x0 < 0) x0 = 0;
				xn = (int) floor((m->x[i] + d) / g->s[0]);
				if (xn >= g->n[0]) xn = g->n[0] - 1;

				k = (l + y) * g->n[0] + x0;
				for (x = x0; x <= xn; ++x) {
					d = x * g->s[0] - m->x[i];
					r2 = dz2dy2 + d * d;

					ra = sqrt(r2);
					q_r = m->q[i] / MAX(ra, RCOR);
					asympcr[k] -= ccutcr * q_r;
					asymphr[k] -= ccuthr * q_r;

					if (ra == 0.0) {
						asympcr[k] += casympcr * m->q[i];
						asymphr[k] += casymphr * m->q[i];
					} else {
						q_r = m->q[i] / ra;
						d = exp(-xappa * ra);
						rs = ra / smear;

						asympcr[k] += q_r * (1.0 - erfc(rs));
						asymphr[k] += 0.5 * q_r * b * (d * erfc(a - rs) - erfc(a + rs) / d);
					}
					++k;
				}
			}
		}
	}
	d = 18.2223 * 18.2223 / w->m.t;
	for (k = 0; k < g->nr; ++k) {
		asympcr[k] *= d;
		asymphr[k] *= d;
	}
}

void asympk(const wvec_t *wvec, const water_t *w, const mol_t *m, double th, double *asympck, double *asymphk, double *huvk0)
{
	int i, k, l;
	double p, sumsin, sumcos, pc, ph, a;
	double smear, xappa, smear2_4, xappa2;
	double *phase, *sinp, *cosp;
	smear = 1.0;
	xappa = 0.0;
	smear2_4 = 0.25 * smear * smear;
	xappa2 = xappa * xappa;
	a = 4.0 * M_PI * 18.2223 * 18.2223 / (w->m.t * wvec->volume);
	
	phase = (double *) calloc(3 * m->natom, sizeof(double));
	sinp = phase + m->natom;
	cosp = sinp + m->natom;

	l = 2;
	i = 3;
	for (k = 1; k < wvec->nk2; ++k) {
		/*l = 3 * k;
		for (i = 0; i < m->natom; ++i)
			phase[i] = g->v[l] * m->x[i] + g->v[l + 1] * m->y[i] + g->v[l + 2] * m->z[i];
		*/
		p = a * exp(-smear2_4 * wvec->v2[k]);
		pc = p / wvec->v2[k];
		ph = M_1_2PI2 * p / (wvec->v2[k] + xappa2);

		if (pc < th && ph < th) {
			sumcos = 1.0;
			sumsin = 0.0;
			i += 3;
		} else {
			cblas_dcopy(m->natom, m->x, 1, phase, 1);
			cblas_dscal(m->natom, wvec->v[i++], phase, 1);
			cblas_daxpy(m->natom, wvec->v[i++], m->y, 1, phase, 1);
			cblas_daxpy(m->natom, wvec->v[i++], m->z, 1, phase, 1);

			vec_sincos(m->natom, phase, sinp, cosp);

			sumcos = cblas_ddot(m->natom, cosp, 1, m->q, 1);
			sumsin = cblas_ddot(m->natom, sinp, 1, m->q, 1);
		}
		asympck[l] = pc * sumcos;
		asymphk[l] = ph * sumcos;
		++l;
		asympck[l] = pc * sumsin;
		asymphk[l] = ph * sumsin;
		++l;
	}
	asympck[0] = 0.0;
	asympck[1] = 0.0;
	asymphk[0] = 0.0;
	asymphk[1] = 0.0;

	cblas_dcopy(m->natom, m->x, 1, phase, 1);
	cblas_dscal(m->natom, wvec->v[0], phase, 1);
	cblas_daxpy(m->natom, wvec->v[1], m->y, 1, phase, 1);
	cblas_daxpy(m->natom, wvec->v[2], m->z, 1, phase, 1);

	vec_sincos(m->natom, phase, sinp, cosp);

	sumcos = cblas_ddot(m->natom, cosp, 1, m->q, 1);
	sumsin = cblas_ddot(m->natom, sinp, 1, m->q, 1);

	l = 0;
	for (i = 0; i < w->natom; ++i) {
		huvk0[l++] = 0.0;
		huvk0[l++] = 0.0;
	}
	free(phase);
}

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i, m, ns = 0;
	double den, dif, dift, ho, hp, w;
	double *c, *d;

	dif = fabs(x - xa[1]);
	c = calloc(2 * n, sizeof(double));
	d = c + n;
	for (i = 0; i < n; i++) {
		if ((dift = fabs(x - xa[i])) < dif) {
			ns = i;
			dif = dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y = ya[ns--];
	for (m = 0; m < (n - 1); m++) {
		for (i = 0; i < (n - m); i++) {
			ho = xa[i] - x;
			hp = xa[i + m] - x;
			w = c[i + 1] - d[i];
			den = ho - hp;
			/* if ((den = ho - hp) == 0.0) nrerror("Error in routine polint"); */
			den = w / den;
			d[i] = hp * den;
			c[i] = ho * den;
		}
		*y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
	}
	free(c);
}
