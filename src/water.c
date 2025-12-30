#define _CRT_SECURE_NO_WARNINGS

#include "water.h"
#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

int waterread(const char *fn, water_t *w)
{
	FILE *f;
	int nb, na, i;
	double *hv, *x;

	f = fopen(fn, "rb");
	if (!f) {
		return 1;
	}
	/* drism */
	if (fread(&w->drism, sizeof(double), 1, f) != 1)
		return 2;
	/* name */
	if (fread(&nb, sizeof(int), 1, f) != 1)
		return 2;
	na = nb < 32 ? nb : 31;
	if (fread(w->name, sizeof(char), na, f) != na)
		return 2;
	w->name[nb] = '\0';
	nb -= na;
	if (nb && fseek(f, nb * sizeof(char), SEEK_CUR))
		return 2;
	/* hoh molecule parameters */
	if (fread(&nb, sizeof(int), 1, f) != 1)
		return 2;
	nb *= sizeof(double);
	na = sizeof(struct WaterMolecule);
	if (nb < na) na = nb;
	if (fread(&w->m, na, 1, f) != 1)
		return 2;
	nb -= na;
	if (nb && fseek(f, nb, SEEK_CUR))
		return 2;
	/* skip g & h*/
	if (fread(&na, sizeof(int), 1, f) != 1)
		return 2;
	if (fread(&nb, sizeof(int), 1, f) != 1)
		return 2;
	if (fseek(f, 2 * (nb * na * sizeof(double) + sizeof(int)), SEEK_CUR))
		return 2;
	/* skip scale */
	if (fread(&nb, sizeof(int), 1, f) != 1)
		return 2;
	if (fseek(f, nb * sizeof(double), SEEK_CUR))
		return 2;
	/* skip name */
	if (fread(&nb, sizeof(int), 1, f) != 1)
		return 2;
	if (fseek(f, nb * sizeof(char), SEEK_CUR))
		return 2;
	/* ngrid */
	if (fread(&w->ngrid, sizeof(int), 1, f) != 1)
		return 2;
	/* dr */
	if (fread(&w->dr, sizeof(double), 1, f) != 1)
		return 2;
	/* skip r */
	if (fseek(f, (w->ngrid - 1) * sizeof(double) + sizeof(int), SEEK_CUR))
		return 2;
	/* dt */
	if (fread(&w->dt, sizeof(double), 1, f) != 1)
		return 2;
	/* skip t, c_fwd & c_bwd */
	if (fseek(f, (w->ngrid + 1) * sizeof(double), SEEK_CUR))
		return 2;
	/* skip strhv */
	if (fread(&na, sizeof(int), 1, f) != 1)
		return 2;
	if (fread(&nb, sizeof(int), 1, f) != 1)
		return 2;
	if (fseek(f, na * nb * sizeof(double), SEEK_CUR))
		return 2;
	/* hv */
	if (fread(&na, sizeof(int), 1, f) != 1)
		return 2;
	if (fread(&w->nfun, sizeof(int), 1, f) != 1)
		return 2;
	nb = w->nfun * w->ngrid;
	hv = (double *) calloc(nb, sizeof(double));
	if (!hv) {
		return 3;
	}
	if (fread(hv, sizeof(double), nb, f) != nb) {
		free(hv);
		return 2;
	}
	if (fread(&w->compc, sizeof(double), 1, f) != 1) {
		free(hv);
		return 2;
	}
	if (fread(&w->comph, sizeof(double), 1, f) != 1) {
		free(hv);
		return 2;
	}
	w->xvv = (double *) malloc(3 * (w->ngrid + 1) * sizeof(double));
	if (!w->xvv) {
		free(hv);
		return 4;
	}
	x = w->xvv;
	x[0] = 2.0 * w->comph;
	x++;
	for (i = 0; i < w->ngrid; ++i)
		x[i] = hv[i] + hv[w->ngrid + i];
	x += w->ngrid;
	x[0] = M_SQRT2 * w->comph;
	x++;
	for (i = 0; i < w->ngrid; ++i)
		x[i] = M_SQRT2 * hv[2 * w->ngrid + i];
	x += w->ngrid;
	x[0] = w->comph;
	memcpy(x + 1, hv + 4 * w->ngrid, w->ngrid * sizeof(double));
	free(hv);

	w->natom = 2;
	w->nfun = 3;
	/* H */
	w->q[0] = w->m.q_h;
	w->s[0] = 0.4;
	w->e[0] = 0.046;

	/* O */
	w->q[1] = -2.0 * w->m.q_h;
	w->s[1] = pow(w->m.A_oo / w->m.B_oo, 1.0 / 6.0);
	w->e[1] = 0.25 * w->m.B_oo * w->m.B_oo / w->m.A_oo;

	/* n */
	w->n[0] = 2;
	w->n[1] = 1;

	fclose(f);
	return 0;
}

void rm_water(water_t *w)
{
	free(w->xvv);
}

int mkxvva(const grid_t *g, const water_t *w, double *xvva)
{
	double *cs, *t;
	int ng, i, j;

	if (g->a[g->na - 1] > (M_PI / w->dr))
		return 1;

	/* dt = M_PI / (w->ngrid * w->dr); */
	ng = w->ngrid + 1;
	t = (double *) calloc(g->na + 4 * (ng - 1), sizeof(double));
	if (!t)
		return 2;

	cs = t + g->na;
	for (i = 0; i < g->na; ++i)
		t[i] = g->a[i] / w->dt;
	/* spline interpolation */
	i = 0;
	for (j = 0; j < 3; ++j) {
		cubicuni(ng, w->xvv + j * ng, cs, FDERIV, 0.0, SDERIV, 0.0);
		resample(ng-1, cs, g->na, t, xvva + j * g->na);
	}
	free(t);

	return 0;
}
