#ifndef __RISM3D_WATER_H
#define __RISM3D_WATER_H

struct WaterMolecule
{
	double q_h;
	double t;
	double rho;
	double eps;
	double r_hh;
	double r_oh;
    double A_oo;
	double B_oo;
	double A_oh;
	double B_oh;
	double A_hh;
	double B_hh;
	double ngc;
};
struct Water
{
	double drism;
	char name[32];
	struct WaterMolecule m;
	int ngrid;
	int natom;
	int nfun;
	double dr, dt;
	double *xvv;
	double compc;
	double comph;
	int n[2];
	double s[2], q[2], e[2];
};
typedef struct Water water_t;

#include "grid.h"

int waterread(const char *, water_t *);
int mkxvva(const grid_t *, const water_t *, double *);

#endif //__RISM3D_WATER_H
