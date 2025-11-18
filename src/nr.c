#include <memory.h>
#include <assert.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>
#include <time.h>
#include <stdlib.h>

#include "eqoz.h"

#define print   printf
#define flush	fflush(stdin)

#define narm    (5)
#define sarm    (0.618033988749895)


/* + 5 * n * sizeof(float)*/
int bicgstab(lneq_t *eq, int N, const float *b, float *x, float *tol, int *it)
{
    float *r, *h, *p, *v, *s;
    float rTh, rTr, normr, norms, alpha, beta, omega, sTr;
    float err, normb;
    int its, maxit;
  
    r = malloc((size_t) N * sizeof(float));
    h = malloc((size_t) N * sizeof(float));
    p = malloc((size_t) N * sizeof(float));
    v = malloc((size_t) N * sizeof(float));
    s = malloc((size_t) N * sizeof(float));
    
    its = 0;
    maxit = *it * 2;
    normb = cblas_snrm2(N,b,1);
    if (normb == 0.0f)
        normb = 1.0f;
    err = *tol * normb;
    
    Jx(eq, x, r);
    ++its;

    cblas_saxpy(N,-1.,b,1,r,1); /* r = Ax-b */
    cblas_scopy(N,r,1,h,1);
    cblas_scopy(N,r,1,p,1);
    
    rTh=cblas_sdot(N,r,1,h,1); /* rho1 */
    normr=cblas_snrm2(N,r,1);
    if (normr <= err) {
        *tol = normr / normb;
        *it = its;
        free(r); free(h); free(p); free(v); free(s);
        return 0;
    }
    while (its < maxit) {
        Jx(eq,p,v);  /* Jx(eq, p, v) */
        ++its;
        
        alpha=rTh/cblas_sdot(N,h,1,v,1);
        cblas_scopy(N,r,1,s,1);  /* s = h */
        cblas_saxpy(N,-alpha,v,1,s,1); /* s = s - alpha * v */
        norms=cblas_snrm2(N,s,1);
        if (norms <= err) {
            cblas_saxpy(N,-alpha,p,1,x,1);
            *tol = norms / normb;
            *it = its;
            free(r); free(h); free(p); free(v); free(s);
            return 0;
        }
        
        Jx(eq,s,r);
        ++its;
        sTr=cblas_sdot(N,s,1,r,1);
        rTr=cblas_sdot(N,r,1,r,1);
        if ( fabs(sTr)<1e-40 || fabs(rTr)<1e-40 )
            omega = 0.;
        else
            omega = sTr/rTr;
        
        cblas_saxpy(N,-alpha,p,1,x,1);
        cblas_saxpy(N,-omega,s,1,x,1);
        cblas_sscal(N,-omega,r,1);
        cblas_saxpy(N,1.f,s,1,r,1);

        normr=cblas_snrm2(N,r,1);
        if (normr <= err) {
            *tol = normr / normb;
            *it = its;
            free(r); free(h); free(p); free(v); free(s);
            return 0;
        }
        if (omega == 0.0f) {
            *tol = normr / normb;
            *it = its;
            free(r); free(h); free(p); free(v); free(s);
            return 2;            
        }
      
        beta=(alpha/omega)/rTh;
        rTh=cblas_sdot(N,r,1,h,1);
        if (rTh == 0.0f) {
            *tol = normr / normb;
            *it = its;
            free(r); free(h); free(p); free(v); free(s);
            return 3;            
        }
        beta*=rTh;
             
        cblas_sscal(N,beta,p,1);
        cblas_saxpy(N,1.,r,1,p,1);
        cblas_saxpy(N,-beta*omega,v,1,p,1);
    }
    free(r); free(h); free(p); free(v); free(s);
  
    *tol = rTr / normb;
    *it = its;
    return 1; /* no convergence */
}

/* Eisenstat, Walker */
typedef struct
{
    double etaM;
    double eta;
    double gamma;
    double prec;
    double normdT;
} fterm_ew_t;

void ftew_setdef(fterm_ew_t *f, double tol)
{
    f->etaM = 0.9;
    f->eta = sqrt(f->etaM / 0.9);
    f->gamma = 0.9;
    f->prec = tol;
    f->normdT = 1.0;
}

double ftew_next(fterm_ew_t *f, double normd)
{
    double etaT;
    etaT = f->gamma * f->eta * f->eta;
    f->eta = normd / f->normdT;
    f->eta = f->gamma * f->eta * f->eta;
    if (etaT > .1 && etaT > f->eta)
        f->eta = etaT;
    if (f->etaM < f->eta)
        f->eta = f->etaM;
    etaT = f->prec / normd;
    if (etaT > f->eta)
        f->eta = etaT;
    f->normdT = normd;
    
    return f->eta;
}

/* 9 * n * sizeof(float) */
int nr(eqoz_t *eq, double *t, double *tol, int *maxit)
{
    double *d, *s, *tN, norms, eta, lambda, err, errN, v;
    float *b, *x, eps;
    int n, m, i, j, flag, nlit, maxlit;
    double lntm, nrtm, lntmi;
    lneq_t ln;
    fterm_ew_t fterm;

    nrtm = clock() / (double) CLOCKS_PER_SEC;
    lntm = .0;
    
    ftew_setdef(&fterm, *tol);
    maxlit = 1000;
    
    ln.natv = eq->solv.natv;
	ln.symc = eq->solv.symc;
    ln.lxvva = eq->solv.lxvva;
    ln.xvva = eq->solv.xvva;
    /*
    m = ln.natv * ln.natv * ln.lxvva;
    ln.xvva = malloc(m * sizeof(float));
    for (j = 0; j < m; ++j)
        ln.xvva[j] = (float) eq->solv.xvva[j];
     */
    ln.indga = eq->solv.indga;
    ln.grid = &eq->lngr;
    m = ln.grid->nr * ln.natv;
    ln.dcdg = (float *) malloc(m * sizeof(float));
    b = (float *) malloc((size_t) m * sizeof(float));
    x = (float *) malloc((size_t) m * sizeof(float));
    
	n = eq->solv.natv * eq->grid.nr;
    v = sqrt(n);
    
    d = (double *) malloc((size_t) n * sizeof(double));
    s = (double *) malloc((size_t) n * sizeof(double));
    tN = (double *) malloc((size_t) n * sizeof(double));
    
    eqoz(&eq->grid, &eq->rism, &eq->solv, eq->closure, t, d, ln.dcdg);
    err = cblas_dnrm2(n, d, 1) / v;
    
    print("||Z(t0)|| = %.1e\n", err);
    print("    #     Nit F      t,s M   eta   ||Jx-Z||       ||Z||\n");
	flush;
    i = 0;
    while (i <= *maxit && err > *tol) {
        
        nlit = maxlit;
        eta = ftew_next(&fterm, err);
        eps = (float) eta;
        for (j = 0; j < m; ++j) {
            b[j] = (float) d[j];
            x[j] = .0f;
        }
        lntmi = clock() / (double) CLOCKS_PER_SEC;
        flag = bicgstab(&ln, n, b, x, &eps, &nlit);
        lntmi = clock() / (double) CLOCKS_PER_SEC - lntmi;
        lntm += lntmi;
        
        if (flag) {
            /* no convergence of linear solver */
			return 1;
        }
        norms = (double) eps;
        for (j = 0; j < m; ++j) {
            s[j] = (double) -x[j];
        }

        
        cblas_dcopy(n, t, 1, tN, 1);
        cblas_daxpy(n, 1.0, s, 1, tN, 1);
        eqoz(&eq->grid, &eq->rism, &eq->solv, eq->closure, tN, d, ln.dcdg);
        errN = cblas_dnrm2(n, d, 1) / v;

        j = 0;
        lambda = 1.0;
        while (j < narm && errN > err) {
            lambda *= sarm;
            cblas_dcopy(n, t, 1, tN, 1);
            cblas_daxpy(n, lambda, s, 1, tN, 1);
            eqoz(&eq->grid, &eq->rism, &eq->solv, eq->closure, tN, d, ln.dcdg);
            errN = cblas_dnrm2(n, d, 1) / v;
            ++j;
        } 
        
        if (j >= narm) {
            /* no convergence of armigo*/
			return 2;
        }
        cblas_dcopy(n, tN, 1, t, 1);
        err = errN;
        
        print("%5d %7d %1d %8.1f %1d %5.2f %10.6f %11.6f\n", i+1, nlit, flag, lntmi, j, eta, norms, err);
		flush;
        ++i;
    }
    
    free(ln.dcdg);
    free(b);
    free(x);
    free(d);
    free(s);
    free(tN);
    
    *maxit = i;
    *tol = err;

    nrtm = clock() / (double) CLOCKS_PER_SEC - nrtm;
    print("NR time is %.1fs, BiCGStab time is %.1fs\n", nrtm, lntm);
    print("||Z(tn)|| = %.1e\n", err);
    
    return 0; /* OK */
}

