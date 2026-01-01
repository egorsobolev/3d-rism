#include "grid.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

void grid_init(const box_t *b, grid_t *g)
{
	int i;
	for (i = 0; i < 3; ++i) {
		/* copy box parameters */
		g->l[i] = b->l[i];
		g->n[i] = b->n[i];
		g->s[i] = b->s[i];
		g->fft_shape[i] = b->n[2 - i];
	}
	g->nr = g->n[0] * g->n[1] * g->n[2];
	g->nk2 = (g->n[0] / 2 + 1) * g->n[1] * g->n[2];
	g->nk = 2 * g->nk2;
	g->volume = g->l[0] * g->l[1] * g->l[2];
}

int mk_wavevectors(grid_t *g, wvec_t *wvec)
{
	int n[] = { g->n[0] / 2 + 1, g->n[1], g->n[2] };
	double a, dk, *k, *kx, *ky, *kz, *g2, *gv, *ga;
	int l, i, j, x, y, z, k0, *ig2, *iga, *index_fwd;
	size_t nk;

	nk = g->nk2;

	/* local arrays
	 kx: double[n0] - wave coord x
	 ky: + double[n1] - wave coord y
	 kz: + double[n2] - wave coord z
	 ig2: + int[nk] - ascending index full / inverse index
	 index_fwd: + int[nk] - forward index (original -> unique)
	 */
	k = (double *) malloc((n[0] + n[1] + n[2]) * sizeof(double) +
	                      2 * nk * sizeof(int));
	if (!k) return -1;
	kx = k;
	ky = kx + n[0];
	kz = ky + n[1];

	ig2 = (int *) (kz + n[2]);
	index_fwd = ig2 + nk;

	/* wave vectors
	 g2: double[nk] (short v2) - squared wave vector
	 gv: + double[3*nk] (short v) - wave vector
	 */
	g2 = (double *) calloc(4 * nk, sizeof(double));
	if (!g2) {;
		free(k);
		return -1;
	}
	gv = g2 + nk;

	/* preprocess */
	l = 0;
	for (j = 0; j < 3; ++j) {
		/* make wave coordinates */
		dk = 2. * M_PI / g->l[j];
		k0 = n[j] / 2 - 1;
		for (i = 0; i < n[j]; ++i) {
			k[l++] = dk * (k0 - (k0 + i) % (int) g->n[j]);
		}
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
	index_fwd[ig2[0]] = 0;
	for (i = 1; i < nk; ++i) {
		if (g2[ig2[i]] != a) {
			++l;
			a = g2[ig2[i]];
			ig2[l] = ig2[i];  // unique -> original
		}
		index_fwd[ig2[i]] = l;  // forward index (original -> unique)
	}
	++l;

	/* make sorted unique array of wave lengths 
	 ga: double[na] - modules of wave vectors
	 iga: int[nk] - index
	 */
	ga = (double *) malloc(l * sizeof(double) + nk * sizeof(int));
	if (!ga) {
		free(k);
		free(g2);
		return -1;
	}
	iga = (int *) (ga + l);

	memcpy(iga, index_fwd, nk * sizeof(int));
	for (i = 0; i < l; i++)
		ga[i] = sqrt(g2[ig2[i]]);

	/* fill Grid structure */
	wvec->nk2 = nk;
	wvec->na = l;
	wvec->volume = g->volume;
	wvec->a = ga;
	wvec->ia = iga;
	wvec->v = gv;
	wvec->v2 = g2;

	/* deallocate memory */
	free(k);

	return 0;
}

void null_wavevectors(wvec_t *wvec)
{
	wvec->a = wvec->v2 = NULL;
}

void rm_wavevectors(wvec_t *wvec)
{
	free(wvec->a);
	free(wvec->v2);
}
