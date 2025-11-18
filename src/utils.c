#include "utils.h"

int cubicuni(int m, double *f, double *cs, int ltype, double lval, int rtype, double rval)
{
	int i, n, n1;
	double df2, df1, c1, c2, d;
	double *A, *B, *C, *D;
	
	if (m < 2)
		return 1;
	n = m - 1;
	n1 = n - 1;
	A = cs;
	B = A + n;
	C = B + n;
	D = C + n;
	/* Si(x) = Ai + Bi(x-xi) + Ci(x-xi)^2 + Di(x-xi)^3
	 * решение 3-х диагональной системы на коэффициенты Ci
	 * первое уравнение - левое граничное условие
	 */

	df1 = f[1] - f[0];
	switch (ltype) {
	case PARABOLA:
		A[0] = 1.0;
		B[0] = 0.0;
		break;
	case FDERIV:
		A[0] = -0.5;
		B[0] = 1.5 * (df1 - lval);
		break;
	case SDERIV:
		A[0] = 0.0;
		B[0] = 0.5 * lval;
		break;
	default:
		return 2;
	}
	/* прямая прогонка перестановки */
	for (i = 1; i < n; ++i) {
		df2 = f[i + 1] - f[i];
		d = 4.0 + A[i - 1];
		A[i] = -1.0 / d;
		B[i] = (3.0 * (df2 - df1) - B[i - 1]) / d;
		df1 = df2;
	}
	/* последнее уравнение - правое граничное условие */
	switch (rtype) {
	case PARABOLA:
		c2 = B[n1] / (1.0 - A[n1]);
		break;
	case FDERIV:
		c2 = (3.0 * (rval - df1) - B[n1]) / (2.0 + A[n1]);
		break;
	case SDERIV:
		c2 = 0.5 * rval;
		break;
	default:
		return 3;
	}
	/* обратная прогонка */
	for (i = n1; i >= 0; --i) {
		c1 = A[i] * c2 + B[i];
		A[i] = f[i];
		B[i] = f[i + 1] - f[i] - (2.0 * c1 + c2) / 3.0;
		C[i] = c1;
		D[i] = (c2 - c1) / 3.0;
		c2 = c1;
	}
	return 0;
}

void resample(int n, const double *cs, int m, const double *x, double *f)
{
	int i, j, n1;
	double t, xi;
	const double *A, *B, *C, *D;

	n1 = n - 1;
	A = cs;
	B = A + n;
	C = B + n;
	D = C + n;

	for (i = 0; i < m; ++i) {
		xi = x[i];
		j = xi < n1 ? (int) xi : n1;
		t = xi - j;
		f[i] = A[j] + t * (B[j] + t * (C[j] + t * D[j]));
	}
}

#define CUTOFF 8
/*#define STKSIZ (8 * sizeof(void*) - 2)*/
#define STKSIZ (8 * sizeof(void*))

void quicksort(int n, double *a, int *idx)
{
	int lostk[STKSIZ];
	int histk[STKSIZ];
	int p, hi, lo, i, j, y, k, l;
	double x;

	if (n < 2)
		return;

	lostk[0] = 0;
	histk[0] = n - 1;
	p = 1;

	while (p) {
		--p;
		lo = lostk[p];
		hi = histk[p];
		i = lo;
		j = hi;
		l = k = (lo + hi + 1) / 2;
		x = a[idx[k]];
		while (i < l || k < j) {
			while (i < l && a[idx[i]] <= x) {
				if (a[idx[i]] == x) {
					--l;
					y = idx[i];
					idx[i] = idx[l];
					idx[l] = y;
				} else
					++i;
			}
			while (k < j && x <= a[idx[j]]) {
				if (x == a[idx[j]]) {
					++k;
					y = idx[j];
					idx[j] = idx[k];
					idx[k] = y;
				} else
					--j;
			}
			if (i < l || k < j) {
				y = idx[i];
				idx[i] = idx[j];
				idx[j] = y;
				if (l == i) ++l;
				if (k == j) --k;
				/*++i;
				--j;*/
			}
		}
		--i;
		if (lo < i) {
			lostk[p] = lo;
			histk[p] = i;
			++p;
		}
		++j;
		if (j < hi) {
			lostk[p] = j;
			histk[p] = hi;
			++p;
		}
	}
}
