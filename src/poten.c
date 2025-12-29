#include "poten.h"

#include <stdlib.h>
#include <memory.h>

#define _USE_MATH_DEFINES
#include <math.h>

#ifndef isnan
#include <float.h>
#define isnan	_isnan
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON	(2.2204460492503131e-016)
#endif

#define RCOR2		(0.050396841995795)
#define RCOR		(0.002)
#define MAX(a,b)	(a > b ? a : b);

void uljuv(const grid_t *g, const water_t *w, const mol_t *m, double rcut, double *u)
{
	int i, j, k, kz, kzy, l, n;
	int x, x0, xn, y, y0, yn, z, z0, zn;
	double *s2, *e, a, d, dz2, dz2dy2, r2, r2s2, s6r6, *ucut, *scut, rcut2;

	rcut2 = rcut * rcut;

	n = w->natom * m->natom;
	ucut = (double *) calloc(2 * n + 2 * w->natom, sizeof(double));
	scut = ucut + w->natom;
	s2 = scut + w->natom;
	e = s2 + n;
	k = 0;
	for (i = 0; i < m->natom; ++i) {
		for (j = 0; j < w->natom; ++j) {
			a = 0.5 * (w->s[j] + m->s[i]);
			s2[k] = a * a;
			e[k] = 4.0 * sqrt(w->e[j] * m->e[i]) / w->m.t;
			++k;
		}
	}

	memset(ucut, 0, 2 * w->natom * sizeof(double));
	l = 0;
	for (i = 0; i < m->natom; ++i) {
		if (rcut >= DBL_EPSILON) {
			for (j = 0; j < w->natom; ++j) {
				r2s2 = rcut2 / s2[l + j];
				s6r6 = 1.0 / (r2s2 * r2s2 * r2s2);
				ucut[j] = e[l + j] * s6r6 * (s6r6 - 1.0);
				scut[j] += ucut[j];
			}
			z0 = (int) ceil((m->z[i] - rcut) / g->s[2]);
			if (z0 < 0) z0 = 0;
			zn = (int) floor((m->z[i] + rcut) / g->s[2]);
			if (zn >= g->n[2]) zn = g->n[2] - 1;
		} else {
			z0 = 0;
			zn = g->n[2] - 1;
		}
		for (z = z0; z <= zn; ++z) {
			d = z * g->s[2] - m->z[i];
			dz2 = d * d;

			if (rcut >= DBL_EPSILON) {
				d = sqrt(rcut2 - dz2);
				y0 = (int) ceil((m->y[i] - d) / g->s[1]);
				if (y0 < 0) y0 = 0;
				yn = (int) floor((m->y[i] + d) / g->s[1]);
				if (yn >= g->n[1]) yn = g->n[1] - 1;
			} else {
				y0 = 0;
				yn = g->n[1] - 1;
			}
			kz = z * g->n[1];
			for (y = y0; y <= yn; ++y) {
				d = y * g->s[1] - m->y[i];
				dz2dy2 = dz2 + d * d;

				if (rcut >= DBL_EPSILON) {
					d = sqrt(rcut2 - dz2dy2);
					x0 = (int) ceil((m->x[i] - d) / g->s[0]);
					if (x0 < 0) x0 = 0;
					xn = (int) floor((m->x[i] + d) / g->s[0]);
					if (xn >= g->n[0]) xn = g->n[0] - 1;
				} else {
					x0 = 0;
					xn = g->n[0] - 1;
				}
				kzy = (kz + y) * g->n[0];
				for (x = x0; x <= xn; ++x) {
					d = x * g->s[0] - m->x[i];
					r2 = dz2dy2 + d * d;

					k = kzy + x;
					for (j = 0; j < w->natom; j++) {
						r2s2 = r2 / s2[l + j];
						if (r2s2 < RCOR2)
							r2s2 = RCOR2;
						s6r6 = 1.0 / (r2s2 * r2s2 * r2s2);

						u[k] += e[l + j] * s6r6 * (s6r6 - 1.0) - ucut[j];
						k += g->nr;
					}
				}
			}
		}
		l += w->natom;
	}
	/*
	for (j = 0; j < w->natom; ++j) {
		printf(" UC%d: %f", j, scut[j]);
	}
	printf("\n");
	*/

	free(ucut);
}

void ucolu(const grid_t *g, const water_t *w, const mol_t *m, const int *spd, double rcut, double *u)
{
	int i, j, k, l;
	int x, x0, xn, y, y0, yn, z, z0, zn;
	int k0, k1, k00, k01, k10, k11;
	double d, dx, dy, dz, dz2, dz2dy2, r2, ra, rcut2;
	double spx, spy, spz, spr[3], sps[3], *spu;
	double x000, x001, x010, x011, x100, x101, x110, x111;
	int spn[3], nsp, spzy;

	rcut2 = rcut * rcut;

	memset(u, 0xff, g->nr * sizeof(double));
	/* prepare sparse grid data structures */
	nsp = 1;
	for (j = 0; j < 3; ++j) {
		spn[j] = g->n[j] / spd[j] + 1;
		sps[j] = spd[j] * g->s[j];
		spr[j] = 1.0 / spd[j];
		nsp *= spn[j];
	}
	spu = (double *) calloc(nsp, sizeof(double));

	/* compute Coloumb potential on sparse grid */
	for (i = 0; i < m->natom; ++i) {
		k = 0;
		for (z = 0; z < spn[2]; ++z) {
			dz = z * sps[2] - m->z[i];
			dz2 = dz * dz;

			for (y = 0; y < spn[1]; ++y) {
				dy = y * sps[1] - m->y[i];
				dz2dy2 = dz2 + dy * dy;

				for (x = 0; x < spn[0]; ++x) {
					dx = x * sps[0] - m->x[i];
					r2 = dz2dy2 + dx * dx;

					ra = sqrt(r2);
					if (ra < RCOR)
						ra = RCOR;

					spu[k] += m->q[i] / ra;
					++k;
				}
			}
		}
	}
	/* copy Coloumb potential from sparse grid to dense grid */
	for (z = 0; z < g->n[2]; z += spd[2]) {
		k0 = g->n[1] * z;
		k1 = spn[1] * (z / spd[2]);
		for (y = 0; y < g->n[1]; y += spd[1]) {
			k00 = g->n[0] * (y + k0);
			k11 = spn[0] * (y / spd[1] + k1);
			for (x = 0; x < g->n[0]; x += spd[0])
				u[x + k00] = spu[x / spd[0] + k11];
		}
	}
	if ((spd[0] * spd[1] * spd[2]) == 1) {
		free(spu);
		return;
	}
	/* compute Coloumb potential within cutoff radius 
	 * for missed points on dense grid
	 */
	for (i = 0; i < m->natom; ++i) {
		r2 = rcut * rcut;

		z0 = (int) ceil((m->z[i] - rcut) / g->s[2]);
		if (z0 < 0) z0 = 0;
		zn = (int) floor((m->z[i] + rcut) / g->s[2]);
		if (zn >= g->n[2]) zn = g->n[2] - 1;

		for (z = z0; z <= zn; ++z) {
			dz = z * g->s[2] - m->z[i];
			dz2 = dz * dz;

			d = sqrt(rcut2 - dz2);
			y0 = (int) ceil((m->y[i] - d) / g->s[1]);
			if (y0 < 0) y0 = 0;
			yn = (int) floor((m->y[i] + d) / g->s[1]);
			if (yn >= g->n[1]) yn = g->n[1] - 1;

			l = z * g->n[1];
			for (y = y0; y <= yn; ++y) {
				dy = y * g->s[1] - m->y[i];
				dz2dy2 = dz2 + dy * dy;

				d = sqrt(rcut2 - dz2dy2);
				x0 = (int) ceil((m->x[i] - d) / g->s[0]);
				if (x0 < 0) x0 = 0;
				xn = (int) floor((m->x[i] + d) / g->s[0]);
				if (xn >= g->n[0]) xn = g->n[0] - 1;

				k = (l + y) * g->n[0];

				spzy = !(z % spd[2]) && !(y % spd[1]);
				x = x0 + (spzy && !(x0 % spd[0]));
				while (x <= xn) {
					if (isnan(u[k + x])) {
						u[k + x] = 0.0;
						for (j = 0; j < m->natom; ++j) {
							dx = x * g->s[0] - m->x[j];
							dy = y * g->s[1] - m->y[j];
							dz = z * g->s[2] - m->z[j];
							r2 = dx * dx + dy * dy + dz * dz;
							ra = sqrt(r2);
							if (ra < RCOR)
								ra = RCOR;

							u[k + x] += m->q[j] / ra;
						}

					}
					++x;
					x += (spzy && !(x % spd[0]));
				}
			}
		}
	}
	/* interpolate */
	for (z = 0; z < g->n[2]; ++z) {
		j = z / spd[2];
		spz = z * spr[2] - j;

		k0 = j * spn[1];
		k1 = k0 + spn[1];

		l = z * g->n[1];
		for (y = 0; y < g->n[1]; ++y) {
			j = y / spd[1];
			spy = y * spr[1] - j;

			k00 = (j + k0) * spn[0];
			k10 = k00 + spn[0];
			k01 = (j + k1) * spn[0];
			k11 = k01 + spn[0];

			k = (l + y) * g->n[0];

			spzy = !(z % spd[2]) && !(y % spd[1]);
			x = spzy && !(x0 % spd[0]);
			while (x < g->n[0]) {
				if (isnan(u[k + x])) {
					j = x / spd[0];
					spx = x * spr[0] - j;

					/* BLEND_103 extends scalar point data into a cube */
					x000 = spu[j + k00];
					x001 = spu[j + k01];
					x010 = spu[j + k10];
					x011 = spu[j + k11];
					++j;
					x100 = spu[j + k00];
					x101 = spu[j + k01];
					x110 = spu[j + k10];
					x111 = spu[j + k11];
					u[k + x] = x000
						+ spx * (x100 - x000 + spy * (x000 - x100 - x010 + x110))
						+ spy * (x010 - x000 + spz * (x000 - x010 - x001 + x011))
						+ spz * (x001 - x000 + spx * (x000 - x100 - x001 + x101))
						+ spx * spy * spz * (-x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111);
				}
				++x;
				x += spzy && !(x % spd[0]);
			}
		}
	}
	free(spu);
}

void ucoluv(const water_t *w, int n, const double *u, double *uuv)
{
	int j, k, l;
	double d;

	d = 18.2223 * 18.2223 / w->m.t;
	l = 0;
	for (j = 0; j < w->natom; ++j) {
		for (k = 0; k < n; ++k) {
			uuv[l + k] += d * w->q[j] * u[k];
		}
		l += n;
	}
}
