#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "jvxl.h"

#ifdef WIN32
#include <crtdefs.h>
#define itoa _itoa
#endif

static unsigned eidx[] =
{
    0, 3, 8, 13, 21, 26, 34, 42, 54
};
static unsigned eset[] =
{
     8,  3,  0,
    11,  8,  3, 2, 0,
     8,  7,  4, 3, 0,
    11,  8,  7, 6, 4, 3, 2, 0,
     9,  8,  3, 1, 0,
    11, 10,  9, 8, 3, 2, 1, 0,
     9,  8,  7, 5, 4, 3, 1, 0,
    11, 10,  9, 8, 7, 6, 5, 4, 3, 2, 1, 0
};
static unsigned A[] = 
{
    0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3
};
static unsigned B[] =
{
    1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7
};
static unsigned ox[] = { 0, 1, 1, 0, 0, 1, 1, 0 };
static unsigned oy[] = { 0, 0, 0, 0, 1, 1, 1, 1 };
static unsigned oz[] = { 0, 0, 1, 1, 0, 0, 1, 1 };

static char range = 90;
static char base = 35;

char *jvxl_edges(const unsigned *n, const unsigned *h, const double *v, double l)
{
	int x, y, z;
    unsigned k, i, j, a, b;
    double va, vb, f;
    char *d, c;
	ptrdiff_t p;
    size_t size;
    
    size = 4096;
    d = (char *) malloc(size);
	if (!d) return NULL;
    p = 0;

    j = 4u;
    for (x = n[2]-2*h[2]; x >= 0; x -= h[2]) {
        j |= 2u;
        for (y = n[1]-2*h[1]; y >= 0; y -= h[1]) {
            j |= 1u;
            for (z = n[0]-2*h[0]; z >= 0; z -= h[0]) {
                for (k = eidx[j]; k < eidx[j + 1]; ++k) {
                    i = eset[k];
                    a = A[i];
                    b = B[i];
                    va = v[x + y + z + ox[a] * h[2] + oy[a] * h[1] + oz[a] * h[0]];
                    vb = v[x + y + z + ox[b] * h[2] + oy[b] * h[1] + oz[b] * h[0]];
                    f = (l - va) / (vb - va);
                    if ((0.0 <= f) && (f < 1.0)) {
                        c = (char) (f * range + base);
                        if (c == 92)
                            c = 33;
                        d[p++] = c;
                        if ((size_t) (p + 1) >= size) {
                            size += 4096;
                            d = realloc(d, size);
							if (!d) return NULL;
						}
                    }
                }
                j &= 6u;                
            }
            j &= 5u;
        }
        j &= 3u;
    }
    d[p] = 0;
    return d;
}

char *jvxl_vertices(const unsigned *n, const unsigned *h, const double *v, double l)
{
    double a, f;
    unsigned s, x, y, z;
    char *d;
	ptrdiff_t p;
    size_t size;
    size = 4096;
    d = malloc(size);
	if (!d) return NULL;

	p = 0;
    f = -1.0;
    s = 0;
    for (x = 0; x < n[2]; x += h[2])
        for (y = 0; y < n[1]; y += h[1])
            for (z = 0; z < n[0]; z += h[0]) {
                a = v[x + y + z] - l;
                if ((f * a) > .0) {
                    ++s;
                } else {
                    f = -f;
                    d[p++] = ' ';
                    p += sprintf(d + p, "%d", s);
                    /*
                    itoa(s, d + p, 10);
                    p += strlen(d + p);
                    */
                    if ((size_t) (p + 12) >= size) {
                        size += 4096;
                        d = realloc(d, size);
						if (!d) return NULL;
                    }
                    s = 1;
                }
            }
    d[p++] = ' ';
    p += sprintf(d + p, "%d", s);
    /*
    itoa(s, d + p, 10);
    p += strlen(d+p);
    */
	return d;
}

int jvxl_write(const char *fn, const jvxl_box_t *b, int nfun, const double *v, const double *lvl)
{
	FILE *f;
	int nc, ne, nv, i, j;
	char *sv, *se;
	unsigned n[3], h[3], t;

    h[2] = 1;
    h[1] = n[2] = (t = b->n[0]);
    h[0] = n[1] = (t *= b->n[1]);
    n[0] = (t *= b->n[2]);

	nc = b->n[2] / 6;
	if (b->n[2] % 6) ++nc;
	nc = (nc + b->n[2] * 13) * b->n[1] * b->n[0];

	f = fopen(fn, "w");
	if (!f)
		return -1;

	fprintf(f, "\n\n%5d%12.6f%12.6f%12.6f [ANGSTROMS]\n", -1, b->o[0], b->o[1], b->o[2]);
	i = 0;
	for (j = 0; j < 3; ++j) {
		fprintf(f, "%5d%12.6f%12.6f%12.6f\n", b->n[j], b->s[i], b->s[i + 1], b->s[i + 2]);
		i += 3;
	}
	fprintf(f, "%5d%12.6f%12.6f%12.6f%12.6f\n", 1, 1.0, 0.0, 0.0, 0.0);
	fprintf(f, "%d 35 90 35 35 Jmol voxel format version 0.9b\n", -nfun);

	j = 0;
	for (i = 0; i < nfun; ++i) {
		sv = jvxl_vertices(n, h, v + j, lvl[i]);
		if (!sv) {
			fclose(f);
			return -2;
		}
		nv = strlen(sv);
		se = jvxl_edges(n, h, v + j, lvl[i]);
		if (!se) {
			fclose(f);
			free(sv);
			return -2;
		}
		ne = strlen(se);
		fprintf(f, "%g %d %d -1 compressionRatio=%g\n%s\n%s\n", lvl[i], nv, ne, (double) nc / (double) (nv + ne), sv, se);
		free(sv);
		free(se);
		j += t;
	}

	fclose(f);
	return 0;
}
