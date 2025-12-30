#ifndef __RISM3D_MOL_H
#define __RISM3D_MOL_H

struct Molecule {
	int natom;
	double *x, *y, *z, *s, *q, *e;
	double o[3];
};
typedef struct Molecule mol_t;

#include "grid.h"

int molread(const char *, mol_t *);
void rm_mol(mol_t *);
void molcenter(const box_t *, mol_t *);

#endif //__RISM3D_MOL_H
