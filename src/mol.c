#include "mol.h"

#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int molread(const char *fn, mol_t *m)
{
	FILE *f;
	char b[20];
	int n, c, i, j, l;
	double *x;
	static int fl[] = {12, 12, 12, 10, 16, 10};

	f = fopen(fn, "rt");
	if (!f) {
		return 1;
	}
	n = 0;
	while (!feof(f) && !n) {
		fgets(b, 9, f);
		l = strlen(b);
		if (b[0] != '#' && l != 0) {
			n = atoi(b);
			if (n < 1)
				return 2;
		}
		if (l == 8) {
			c = fgetc(f);
			/* skip rest of line */
			while (!feof(f) && c != '\n')
				c = fgetc(f);
		}
	}
	x = (double *) malloc(6 * n * sizeof(double));
	if (!x)
		return 3;

	i = 0;
	while (!feof(f) && i < n) {
		c = fgetc(f);
		if (c != '#' && c != '\n') {
			ungetc(c, f);
			for (j = 0; j < 6; ++j) {
				fgets(b, fl[j] + 1, f);
				l = strlen(b);
				if (l < fl[j]) {
					free(x);
					return 2;
				}
				x[j * n + i] = atof(b);
			}
			c = (char) fgetc(f);
			++i;
		}
		/* skip rest of line */
		while (!feof(f) && c != '\n')
			c = (char) fgetc(f);
	}
	fclose(f);
	if (i < n) {
		free(x);
		return 2;
	}
	
	m->natom = n;
	m->x = x;
	m->y = m->x + n;
	m->z = m->y + n;
	m->s = m->z + n;
	m->q = m->s + n;
	m->e = m->q + n;

	return 0;
}

void null_mol(mol_t *m)
{
	m->x = NULL;
}

void rm_mol(mol_t *m)
{
	free(m->x);
}

void molcenter(const box_t *b, mol_t *m)
{
	int i, j;
	double c, *x;

	x = m->x;
	for (j = 0; j < 3; ++j) {
		c = x[0];
		for (i = 1; i < m->natom; ++i) {
			c += x[i];
		}
		c = 0.5 * b->l[j] - c / m->natom;
		m->o[j] = -c;
		for (i = 0; i < m->natom; ++i)
			x[i] += c;

		x += m->natom;
	}
}
