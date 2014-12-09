/******************************************************************************
 *
 * File:           qp.c
 *
 * Created:        07/03/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Quadratic Programming problem solver.
 *
 * Description:    A wrapper to ql0001().
 *
 * Revisions:      None
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
/*#include <values.h>*/
#include <assert.h>
#include <f2c.h>
#include <math.h>
#include <errno.h>
#include <limits.h>

#if !defined(UROUND)
#define UROUND 1.11e-16
#endif
#define EPS 1.0e-10

static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "error: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(1);
}

/* Allocates n1xn2 matrix. Note that it will be accesses as [n2][n1] */
static double** alloc2d(int n1, int n2)
{
    unsigned int size = (unsigned int) n1 * (unsigned int) n2;
    double* p;
    double** pp;
    int i;

    assert(n1 > 0);
    assert(n2 > 0);
    assert((long) n1 * (long) n2 <= UINT_MAX);
    if ((p = calloc(size, sizeof(double))) == NULL)
        quit("dalloc2d(): %s\n", strerror(errno));

    size = n2 * sizeof(double*);
    if ((pp = malloc(size)) == NULL)
        quit("dalloc2d(): %s\n", strerror(errno));
    for (i = 0; i < n2; i++)
        pp[i] = p + i * n1;

    return pp;
}

static void free2d(double** pp)
{
    double* p;

    assert(pp != NULL);
    p = pp[0];
    free((void*) pp);
    assert(p != NULL);
    free((void*) p);
}

/* C -> Fortran matrix convertion */
static void convert(double** Ain, int N, int M, double* Aout)
{
    double* A = Aout;
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++)
            A[j] = Ain[j][i];
        A = &A[M];
    }
}

/* Interface to ql0001_() */
int ql0001(int nc, int nec, int nv, double** C, double* D, double** A, double* B, double* Xmin, double* Xmax, double* X)
{
    /*
     * We treat ql0001_() as a black box here. Since it takes addresses
     * rather than values, we create a variable for each of the input
     * parameters in case it is modified 
     */
    integer M = nc;
    integer MMAX = nc;
    integer ME = nec;
    integer N = nv;
    integer NMAX = N;
    integer MNN = M + N + N;
    integer LWAR = 3 * N * N / 2 + 10 * N + 2 * M + 1;
    integer LIWAR = N;
    integer iout = 0;
    integer ifail = 0;
    integer iprint = 1;
    doublereal eps = UROUND;
    doublereal* U = calloc((nc + nv * 2), sizeof(double));
    doublereal* WAR = calloc(LWAR, sizeof(double));
    integer* IWAR = calloc(N, sizeof(int));
    double* CF = calloc(nv * nv, sizeof(double));
    double* AF = calloc(nv * nc, sizeof(double));

    extern int ql0001_(integer* m, integer* me, integer* mmax, integer* n, integer* nmax, integer* mnn, doublereal* c, doublereal* d, doublereal* a, doublereal* b, doublereal* xl, doublereal* xu, doublereal* x, doublereal* u, integer* iout, integer* ifail, integer* iprint, doublereal* war, integer* lwar, integer* iwar, integer* liwar, doublereal* eps);

    assert(sizeof(int) == sizeof(integer));
    assert(sizeof(double) == sizeof(doublereal));

    /*
     * ! 
     */
    IWAR[0] = 1;

    convert(C, nv, nv, CF);
    convert(A, nv, nc, AF);

    ql0001_(&M, &ME, &MMAX, &N, &NMAX, &MNN, CF, D, AF, B, Xmin, Xmax, X, U, &iout, &ifail, &iprint, WAR, &LWAR, IWAR, &LIWAR, &eps);

    free(AF);
    free(CF);
    free(IWAR);
    free(WAR);
    free(U);

    return (int) ifail;
}

/* libf2c.a wants it */
int MAIN__()
{
    return 0;
}

/** Solver for a Quadratic Programming problem
 * @param nec Number of equality contraints
 * @param nic Number of inequality constraints
 * @param nwc Number of "soft" constraints
 * @param nv Number of variables
 * @param A Data matrix [nec+nic][nv] for equality and inequality constraints
 * @param B Constant vector [nec+nic] for equality and inequality constraints
 * @param Aw Data matrix [nwc][nv] for soft constraints
 * @param Bw Constant vector [nwc] for soft constraints
 * @param W Vector [nwc] of weights for soft constraints in the objective 
 *          function
 * @param Xmin Vector [nv] of lower bound values for the variables
 * @param Xmax Vector [nv] of upper bound values for the variables
 * @param X Vector [nv] of variables; on return contains solution
 * @return 1 on success; 0 on failure
 *
 * A[i] X + B[i] = 0, 0 <= i <= nec - 1
 * A[i] X + B[i] >= 0, nec <= i <= nec + nic - 1
 * <sum, 0 <= i <= nwc - 1> W[i] (Aw[i] X + Bw[i])^2 = min
 * Xmin[i] <= X[i] <= Xmax[i], 0 <= i <= nv - 1
 */
int qp_solve(int nec, int nic, int nwc, int nv, double** A, double* B, double** Aw, double* Bw, double* W, double* Xmin, double* Xmax, double* X)
{
    int nc = nec + nic;
    double** C = alloc2d(nv, nv);
    double* D = calloc(nv, sizeof(double));
    int ok = 1;
    int i, j, c;

    if (nv <= 0)
        quit("qp_solve(): nv <= 0\n");
    if (nwc <= 0)
        quit("qp_solve(): nwc <= 0\n");

    /*
     * check that the bilinear form is positively defined -- this is
     * necessary to obtain a positively defined matrix C in QL0001() -- see
     * qld.f 
     */
    for (c = 0; c < nwc; ++c)
        if (W[c] < 0.0)
            quit("qp_solve(): W[%d] = %.3g < 0\n", c, W[c]);

    /*
     * calculate contribution of each weak constraint into objective
     * function 
     */
    for (c = 0; c < nwc; ++c) {
        double* row = Aw[c];
        double w = W[c];
        double b = Bw[c];

        for (i = 0; i < nv; ++i)
            D[i] += w * b * row[i];

        for (i = 0; i < nv; ++i) {
            double w_row_i = w * row[i];

            if (w_row_i == 0.0)
                continue;
            for (j = 0; j < nv; ++j)
                C[i][j] += w_row_i * row[j];
        }
    }

    ql0001(nc, nec, nv, C, D, A, B, Xmin, Xmax, X);

    /*
     * check solution consistency 
     */

    /*
     * check equality constraints 
     */
    for (c = 0; c < nec; ++c) {
        double sum = 0;

        for (i = 0; i < nv; ++i)
            sum += A[c][i] * X[i];
        if (fabs(sum + B[c]) > EPS) {
            ok = 0;
            goto end;
        }
    }
    /*
     * check inequality constraints 
     */
    for (c = nec; c < nc; ++c) {
        double sum = 0.0;

        for (i = 0; i < nv; ++i)
            sum += A[c][i] * X[i];
        if (sum + EPS + B[c] < 0.0) {
            ok = 0;
            goto end;
        }
    }
    /*
     * check bound values for variables 
     */
    for (i = 0; i < nv; ++i) {
        if (X[i] + EPS < Xmin[i] || X[i] - EPS > Xmax[i]) {
            ok = 0;
            goto end;
        }
    }

  end:
    free(D);
    free2d(C);

    return ok;
}

#undef EPS
#undef UROUND
