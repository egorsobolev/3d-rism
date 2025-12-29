#include "grid.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

void mkbox(const double *gsp, const double *margin, const mol_t *m, box_t *b)
{
	double mn, mx, l, *x;
	int i, j, n, k;

	x = m->x;
	for (j = 0; j < 3; ++j) {
		mn = mx = x[0];
		for (i = 1; i < m->natom; ++i) {
			if (x[i] < mn) mn = x[i];
			if (x[i] > mx) mx = x[i];
		}
		l = (mx - mn) + 2.0 * margin[j];

		k = (int) (l / gsp[j] + 0.5);
		k = k + (k & 1);

		n = 1;
		while (k > 1) {
			if (!(k % 2)) {
				k /= 2;
				n *= 2;
			} else if (!(k % 3)) {
				k /= 3;
				n *= 3;
			} else if (!(k % 5)) {
				k /= 5;
				n *= 5;
			} else {
				k = k * n + 2;
				n = 1;
			}
		}
		b->n[j] = n;
		b->l[j] = n * gsp[j];
		b->s[j] = gsp[j];
		x += m->natom;
	}
}

int ginit(const box_t *b, grid_t *g)
{
	int n[] = { b->n[0] / 2 + 1, b->n[1], b->n[2] };
	double a, dk, *k, *kx, *ky, *kz, *g2, *gv, *ga;
	int l, i, j, nk, x, y, z, k0, *ig2, *iga;

	/* (nx / 2 + 1) * ny * nz */
	/* 2 * (nx / 2 + 1) * ny * nz */
	nk = n[0] * n[1] * n[2];

	k = (double *) calloc(n[0] + n[1] + n[2], sizeof(double));
	if (!k) return -1;
	kx = k;
	ky = kx + n[0];
	kz = ky + n[1];

	ig2 = (int *) calloc(nk, sizeof(int));
	if (!ig2) {
		free(k);
		return -1;
	}

	g2 = (double *) calloc(4 * nk, sizeof(double));
	if (!g2) {
		free(k);
		free(ig2);
		return -1;
	}
	gv = g2 + nk;

	iga = (int *) calloc(nk, sizeof(int));
	if (!iga) {
		free(k);
		free(ig2);
		free(g2);
		return -1;
	}

	/* preprocess */
	l = 0;
	for (j = 0; j < 3; ++j) {
		/* copy box parameters */
		g->l[j] = b->l[j];
		g->n[j] = b->n[j];
		g->s[j] = b->s[j];
		/* make wave coordinates */
		dk = 2. * M_PI / b->l[j];
		k0 = n[j] / 2 - 1;
		for (i = 0; i < n[j]; ++i)
			k[l++] = dk * (k0 - (k0 + i) % b->n[j]);
	}

	/* fill grid by wave vectors */
	l = 0;
	i = 0;
	for (z = 0; z < n[2]; ++z)
		for (y = 0; y < n[1]; ++y)
			for (x = 0; x < n[0]; ++x) {
				gv[l++] = kx[x];
				gv[l++] = ky[y];
				gv[l++] = kz[z];
				g2[i++] = kx[x] * kx[x] + ky[y] * ky[y] + kz[z] * kz[z];
			}

	/* make ascending index of g2 by quicksort algorithm */
	for (i = 0; i < nk; ++i)
		ig2[i] = i;
	quicksort(nk, g2, ig2);
	/* make index regenerating g2 from sorted unique wave lengths */
	l = 0;
	a = g2[ig2[0]];
	for (i = 1; i < nk; ++i) {
		if (g2[ig2[i]] != a) {
			++l;
			a = g2[ig2[i]];
		}
		iga[ig2[i]] = l;
	}
	++l;
	/* make sorted unique array of wave lengths */
	ga = (double *) calloc(l, sizeof(double));
	if (!ga) {
		free(k);
		free(ig2);
		free(g2);
		free(iga);
		return -1;
	}
	for (i = 0; i < nk; ++i)
		ga[iga[i]] = sqrt(g2[i]);

	/* fill Grid structure */
	g->nr = b->n[0] * b->n[1] * b->n[2];
	g->nk2 = nk;
	g->nk = 2 * nk;
	g->na = l;
	g->a = ga;
	g->ia = iga;
	g->v = gv;
	g->v2 = g2;

	/* deallocate memory */
	free(k);
	free(ig2);

	return 0;
}

int mkdgrid(int n, int *nn, dgrid_t *p)
{
	int i;

	p->n = n;
	p->nn = calloc(n, sizeof(size_t));
	if (!p->nn)
		return -1;
	for (i = 0; i < n; ++i)
		p->nn[i] = nn[n - i - 1];

	p->b = p->nn[p->n - 1];
	p->nr = p->b;
	p->a = 2 * (p->b / 2 + 1);
	p->nk = p->a;
	for (i = p->n - 2; i >= 0; --i) {
		p->nr *= p->nn[i];
		p->nk *= p->nn[i];
	}
	p->m = p->nk / p->a;

	p->data = (double*) malloc(2 * sizeof(double) * p->nk);
	if (!p->data) {
		free(p->nn);
		return -2;
	}
	p->tmp = p->data + p->nk;
	return 0;
}

void rmdgrid(dgrid_t *p)
{
	free(p->data);
	free(p->nn);
}

void dgrid_fwd(dgrid_t *p)
{
	fft_r2c(p->n, p->nn, p->data, (complex_double *) p->tmp);
	memcpy(p->data, p->tmp, sizeof(double) * p->nk);
}

void dgrid_bwd(dgrid_t *p)
{
	fft_c2r(p->n, p->nn, (complex_double *) p->data, p->tmp);
	memcpy(p->data, p->tmp, sizeof(double) * p->nr);
}


int mkfgrid(int n, int *nn, fgrid_t *p)
{
	int i;

	p->n = n;
	p->nn = calloc(n, sizeof(int));
	if (!p->nn)
		return -1;
	for (i = 0; i < n; ++i)
		p->nn[i] = nn[n - i - 1];

	p->b = p->nn[p->n - 1];
	p->nr = p->b;
	p->a = 2 * (p->b / 2 + 1);
	p->nk = p->a;
	for (i = p->n - 2; i >= 0; --i) {
		p->nr *= p->nn[i];
		p->nk *= p->nn[i];
	}
	p->m = p->nk / p->a;

	p->data = (float*) malloc(2 * sizeof(float) * p->nk);
	if (!p->data) {
		free(p->nn);
		return -3;
	}
	p->tmp = p->data + p->nk;
	return 0;
}

void rmfgrid(fgrid_t *p)
{
	free(p->data);
	free(p->nn);
}

void fgrid_fwd(fgrid_t *p)
{
	fftf_r2c(p->n, p->nn, p->data, (complex_float *) p->tmp);
	memcpy(p->data, p->tmp, sizeof(float) * p->nk);
}

void fgrid_bwd(fgrid_t *p)
{
	fftf_c2r(p->n, p->nn, (complex_float *) p->data, p->tmp);
	memcpy(p->data, p->tmp, sizeof(float) * p->nr);
}
