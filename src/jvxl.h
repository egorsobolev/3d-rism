#ifndef __RISM3D_JVXL_H
#define __RISM3D_JVXL_H

struct JVXLBOX
{
	double s[9];
	double o[3];
	int n[3];
};
typedef struct JVXLBOX jvxl_box_t;

#include "grid.h"
#include "mol.h"

int jvxl_write(const char *fn, const jvxl_box_t *b, int n, const double *v, const double *lvl);

#endif //__RISM3D_JVXL_H
