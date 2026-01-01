#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <memory.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cblas.h>

#include "mol.h"
#include "water.h"
#include "grid.h"
#include "utils.h"
#include "poten.h"
#include "asymp.h"
#include "eqoz.h"
#include "tdyn.h"
#include "jvxl.h"

#include <argtable2.h>

#ifdef WIN32
#define strcasecmp _stricmp
#endif

int main(int argc, char **argv)
{
	struct arg_file *wfile = arg_file1("w", NULL, NULL, "water file");
	struct arg_str *closure = arg_str0("c", NULL, NULL, "closure: HNC, PLHNC (default: PLHNC)");
	struct arg_dbl *prec = arg_dbl0("t", NULL, NULL, "tolerance (default: 1e-5)");
	struct arg_int *mxit = arg_int0(NULL, "maxi", NULL, "maximum number of iterations (default: 300)");
	struct arg_dbl *step = arg_dbl0("s", NULL, NULL, "common spatial step, A (default: 0.5)");
	struct arg_dbl *spx = arg_dbl0(NULL, "sx", NULL, "spatial step by x, A");
	struct arg_dbl *spy = arg_dbl0(NULL, "sy", NULL, "spatial step by y, A");
	struct arg_dbl *spz = arg_dbl0(NULL, "sz", NULL, "spatial step by z, A");
	struct arg_dbl *margin = arg_dbl0("m", NULL, NULL, "common margin, A (default: 15)");
	struct arg_dbl *mgx = arg_dbl0(NULL, "mx", NULL, "margin by x, A");
	struct arg_dbl *mgy = arg_dbl0(NULL, "my", NULL, "margin by y, A");
	struct arg_dbl *mgz = arg_dbl0(NULL, "mz", NULL, "margin by z, A");
	struct arg_dbl *ljcut = arg_dbl0(NULL, "ljcut", NULL, "LJ potential cutoff radius, A (default: 10)");
	struct arg_dbl *ccut = arg_dbl0(NULL, "ccut", NULL, "Coulomb potential cutoff radius, A (default: 5)");
	struct arg_int *thinn = arg_int0("i", NULL, NULL, "common parameter for thinning Coulomb grid (default: 2)");
	struct arg_int *thx = arg_int0(NULL, "thx", NULL, "parameter for thinning Coulomb grid by x");
	struct arg_int *thy = arg_int0(NULL, "thy", NULL, "parameter for thinning Coulomb grid by y");
	struct arg_int *thz = arg_int0(NULL, "thz", NULL, "parameter for thinning Coulomb grid by z");
	struct arg_dbl *ath = arg_dbl0(NULL, "ath", NULL, "parameter for coarsening asymptotics (default: 1e-7)");
	/*struct arg_file *pfile = arg_file0("p", NULL, NULL, "file with/for prepared equation data");*/
	struct arg_int *iso	= arg_intn("l", NULL, NULL, 0, 10, "isosurface level in percent of maximum (up to 10 levels)");
	struct arg_lit *help  = arg_lit0(NULL, "help", "print this help and exit");
	struct arg_lit *vers  = arg_lit0(NULL, "version", "print version information and exit");
	struct arg_file *rtxt = arg_file1(NULL, NULL, "<molecule>", NULL);
	struct arg_end *end = arg_end(20);
	void *argtable[] = {
		wfile, closure,
		prec, mxit,
		step, spx, spy, spz,
		margin, mgx, mgy, mgz,
		ljcut, ccut,
		thinn, thx, thy, thz,
		ath, /*pfile, */iso,
		help, vers,
		rtxt,
		end
	};
	const char* progname = "rism3d";
	int exitcode=0;
	int nerrors; 
	mol_t m;
	water_t w;
	double gsp[3], mrg[3];
	int spd[3];
	box_t b;
	jvxl_box_t bj;
	int err, i, j, k;
	eqoz_t eq;
	double tol, *tuv, *cuv, td[10], *mx, *lvl, th;
	FILE *f;
	musc_t *musc;
	closure_c_t *closure_c;
	char prefix[20], filename[40];

	int nit, flag, nprefix;

	cuv = tuv = NULL;
	null_water(&w);
	null_mol(&m);
	null_eqoz(&eq);

	if (arg_nullcheck(argtable) != 0)
	{
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n", progname);
		exit(1);
	}
	/* set any command line default values prior to parsing */
	closure->sval[0] = "PLHNC";
	prec->dval[0] = 1e-5;
	mxit->ival[0] = 300;
	step->dval[0] = 0.5;
	margin->dval[0] = 15.0;
	ljcut->dval[0] = 10.0;
	ccut->dval[0] = 10.0;
	thinn->ival[0] = 2;
	ath->dval[0] = 1e-7;

	/* Parse the command line as defined by argtable[] */
	nerrors = arg_parse(argc, argv, argtable);

	if (help->count > 0) {
		printf("Usage: %s", progname);
		arg_print_syntax(stdout, argtable,"\n");
		printf("Solve 3D-RISM equations.\n\n");
		arg_print_glossary(stdout, argtable,"  %-20s %s\n");
		exitcode = 0;
		goto exit1;
	}
	/* special case: '--version' takes precedence error reporting */
	if (vers->count > 0){
		printf("March 2011, Egor Sobolev\n");
		exitcode = 0;
		goto exit1;
	}
	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0) {
		/* Display the error details contained in the arg_end struct.*/
		arg_print_errors(stdout, end, progname);
		printf("Try '%s --help' for more information.\n", progname);
		exitcode = 1;
		goto exit1;
	}

	nprefix = (int) (rtxt->extension[0] - rtxt->basename[0]);
	if (nprefix > 19) nprefix = 19;
	strncpy(prefix, rtxt->basename[0], nprefix);
	prefix[nprefix] = '\0';

	gsp[0] = gsp[1] = gsp[2] = step->dval[0];
	if (spx->count) gsp[0] = spx->dval[0];
	if (spy->count) gsp[1] = spy->dval[0];
	if (spz->count) gsp[2] = spz->dval[0];

	mrg[0] = mrg[1] = mrg[2] = margin->dval[0];
	if (mgx->count) mrg[0] = mgx->dval[0];
	if (mgy->count) mrg[1] = mgy->dval[0];
	if (mgz->count) mrg[2] = mgz->dval[0];

	spd[0] = spd[1] = spd[2] = thinn->ival[0];
	if (thx->count) spd[0] = thx->ival[0];
	if (thy->count) spd[1] = thy->ival[0];
	if (thz->count) spd[2] = thz->ival[0];

	th = ath->dval[0];

	if ((spd[0] * spd[1] * spd[2]) < 1) {
		printf("%s: invalid values of parameters for thinning Coloumb grid.\n", progname);
		exitcode = 1;
		goto exit1;
	}

	if (!strcasecmp((char *) closure->sval[0], "HNC")) {
		eq.closure = &hnc;
		closure_c = &hnc_c;
		musc = &musc_hnc;
	} else if (!strcasecmp((char *) closure->sval[0], "PLHNC")) {
		eq.closure = &plhnc;
		closure_c = &plhnc_c;
		musc = &musc_plhnc;
	} else {
		printf("Unknown closure.\n");
		printf("Try '%s --help' for more information.\n", progname);
		exitcode = 0;
		goto exit1;
	}

	if (molread(rtxt->filename[0], &m)) {
		printf("Cannot read molecule data from file '%s'\n", rtxt->filename[0]);
		exitcode = 2;
		goto exit1;
	}

	if (waterread(wfile->filename[0], &w)) {
		printf("Cannot read water data from file '%s'\n", wfile->filename[0]);
		exitcode = 3;
		goto exit1;
	}
	/* Water force constants from Kovalenko */
	/*
	w.s[0] = 2. * 0.69387933 / pow(2.0, 1.0/6.0);
	w.e[0] = 0.0152;
	w.s[1] = 2. * 1.7683 / pow(2.0, 1.0/6.0);
	w.e[1] = 0.152;
	w.m.t = 298.0 * (1.380658e-23 * 6.0221367e+23) / 4184;
	*/

	mkbox(gsp, mrg, &m, &b);
	printf(" BOX                X        Y        Z\n");
	printf("  L:         %8.2lf %8.2lf %8.2lf\n", b.l[0], b.l[1], b.l[2]);
	printf("  N:         %8d %8d %8d\n", b.n[0], b.n[1], b.n[2]);
	printf("\n");

	molcenter(&b, &m);

	// TODO: check exit code and raise exception
	if (mk_eqoz(&eq, &b, &w, &m, ljcut->dval[0], ccut->dval[0], spd, th)) {
		exitcode = 5;
		goto exit1;
	}

	/*
	if (pfile->count) {
		f = fopen(pfile->filename[0], "wb");
		fwrite(&w.natom, sizeof(int), 1, f);
		fwrite(&g.n, sizeof(int), 3, f);
		fwrite(&g.nr, sizeof(int), 1, f);
		fwrite(&g.nk, sizeof(int), 1, f);
		fwrite(eq.rism.uuv, sizeof(double), w.natom * g.nr + 2 * (g.nr + g.nk + w.natom), f);
		fclose(f);
	}
	*/
	tuv = (double *) calloc(eq.grid.nr * w.natom, sizeof(double));
	if (!tuv) {
		printf("%s: insufficient memory\n", progname);
		exitcode = 4;
		goto exit1;
	}
	i = 0;
	for (j = 0; j < w.natom; ++j)
		for (k = 0; k < eq.grid.nr; ++k)
			tuv[i++] = eq.rism.asymcr[k] * eq.solv.charge[j];

	tol = prec->dval[0];
	nit = mxit->ival[0];

	flag = nr(&eq, tuv, &tol, &nit);
	if (flag) {
		if (flag == 1)
			printf("NR: no convergence of linear solver.\n");
		else if (flag == 2)
			printf("NR: no convergence of armigo.\n");

		goto exit1;
	}
	printf("\n");

	cuv = (double *) calloc(eq.grid.nr * w.natom, sizeof(double));
	if (!cuv) {
		printf("%s: insufficient memory\n", progname);
		exitcode = 1;
		goto exit1;
	}
	closure_c(eq.grid.nr * w.natom, eq.rism.uuv, tuv, cuv);

	cblas_daxpy(eq.grid.nr * w.natom, 1.0, cuv, 1, tuv, 1);

	strcpy(filename, prefix);
	strcpy(filename + nprefix, "_huv.rbin");
	f = fopen(filename, "wb");
	fwrite(&eq.grid.nr, sizeof(int), 1, f);
	fwrite(&w.natom, sizeof(int), 1, f);
	fwrite(tuv, sizeof(double), eq.grid.nr * w.natom, f);
	fclose(f);
	strcpy(filename, prefix);
	strcpy(filename + nprefix, "_cuv.rbin");
	f = fopen(filename, "wb");
	fwrite(&eq.grid.nr, sizeof(int), 1, f);
	fwrite(&w.natom, sizeof(int), 1, f);
	fwrite(cuv, sizeof(double), eq.grid.nr * w.natom, f);
	fclose(f);

	td[0] = etot(&eq.grid, &w, eq.rism.uuv, tuv, cuv);
	td[1] = (*musc)(&eq.grid, &w, eq.rism.uuv, tuv, cuv);
	td[2] = mugf(&eq.grid, &w, eq.rism.uuv, tuv, cuv);
	td[3] = ichi(&eq.grid, &w, eq.rism.uuv, tuv, cuv);
	printf("xene = %.3lf musc = %.3lf mugf = %.3lf pmv = %.3lf\n", td[0], td[1], td[2], td[3]);

	free(cuv);

	if (iso->count) {
		strcpy(filename, prefix);
		strcpy(filename + nprefix, "_guv_");
		lvl = (double *) malloc(2 * w.natom * sizeof(double));
		if (!lvl) {
			printf("%s: insufficient memory\n", progname);
			exitcode = 4;
			goto exit1;
		}
		mx = lvl + w.natom;
		memset(bj.s, 0, 9 * sizeof(double));
		for (j = 0; j < 3; ++j) {
			bj.o[j] = m.o[j];
			bj.n[j] = b.n[j];
			bj.s[4 * j] = b.s[j];
		}
		k = 0;
		for (j = 0; j < w.natom; ++j) {
			tuv[k] += 1.0;
			mx[j] = tuv[k];
			++k;
			for (i = 1; i < eq.grid.nr; ++i) {
				tuv[k] += 1.0;
				if (tuv[k] > mx[j])
					mx[j] = tuv[k];
				++k;
			}
		}
		for (i = 0; i < iso->count; ++i) {
			if (0 <= iso->ival[i] && iso->ival[i] < 100) {
				sprintf(filename + nprefix + 5, "%d", iso->ival[i]);
				/*
				itoa(iso->ival[i], filename + nprefix + 5, 10);
				*/
				strcpy(filename + strlen(filename), ".jvxl");
				for (j = 0; j < w.natom; ++j)
					lvl[j] = 0.01 * mx[j] * (double) iso->ival[i];
				jvxl_write(filename, &bj, w.natom, tuv, lvl);
			} else {
				printf("%s, warning: options set invalid value of isosurfce level %d and will be ignored.\n", progname, iso->ival[i]);
			}
		}
		free(lvl);
	}
exit1:
	free(tuv);
	rm_eqoz(&eq);
	rm_water(&w);
	rm_mol(&m);
	arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0])); 
	exit(exitcode);
}
