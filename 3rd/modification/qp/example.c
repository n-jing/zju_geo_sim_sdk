/******************************************************************************
 *
 * File:           qp.h
 *
 * Created:        12/06/2001
 *
 * Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *
 * Purpose:        Examples of using `qp_solve()'
 *
 *****************************************************************************/

#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "qp.h"

static void example1()
{
    double A0[] = { -1.0, -1.0 };
    double Aw0[] = { 1.0, 0.0, 0.0, 1.0 };
    double* A[1];
    double B[] = { 5.0 };
    double W[] = { 1.0, 4.0 };
    double* Aw[2];
    double Bw[] = { -4.0, -2.0 };
    double Xmin[] = { 0.0, 0.0 };
    double Xmax[] = { 3.0, DBL_MAX };
    double X[] = { 1.5, 1.5 };

    /*
     * initialise 2D matrices 
     */
    A[0] = &A0[0];
    Aw[0] = &Aw0[0];
    Aw[1] = &Aw0[2];

/* *INDENT-OFF* */
    printf("\nQP problem:\n"
	   "  -x[0] - x[1] + 5 >= 0\n"
	   "  (x[0] - 4)^2 + 4(x[1] - 2)^2 = min\n"
	   "  0 <= x[0] <= 3\n"
	   "  0 <= x[1] <= inf\n"
	   "  Analytical solution: X = (3, 2)\n"
	   "  Numerical solution: ");
/* *INDENT-ON* */

    qp_solve(0, 1, 2, 2, A, B, Aw, Bw, W, Xmin, Xmax, X);
    printf("X = (%f, %f)\n", X[0], X[1]);
}

static void example2()
{
    double A0[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, -2.0, -2.0 };
    double Aw0[] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0 };
    double* A[2];
    double B[] = { -5.0, 3.0 };
    double W[] = { 1.0, 1.0, 1.0 };
    double* Aw[3];
    double Bw[] = { -1.0, 0.0, 0.0 };
    double Xmin[] = { -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX };
    double Xmax[] = { DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX };
    double X[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

    A[0] = &A0[0];
    A[1] = &A0[5];
    Aw[0] = &Aw0[0];
    Aw[1] = &Aw0[5];
    Aw[2] = &Aw0[10];

/* *INDENT-OFF* */
    printf("\nQP problem:\n"
	   "  x[0] + x[1] + x[2] + x[3] + x[4] - 5 = 0\n"
	   "  x[2] - 2x[3] - 2x[4] + 3 = 0\n"
	   "  (x[0] - 1)^2 + (x[1] - x[2])^2 + (x[3] - x[4])^2 = min\n"
	   "  Analytical solution: X = (1, 1, 1, 1, 1)\n"
	   "  Numerical solution: ");
/* *INDENT-ON* */

    qp_solve(2, 0, 3, 5, A, B, Aw, Bw, W, Xmin, Xmax, X);
    printf("X = (%f, %f, %f, %f, %f)\n", X[0], X[1], X[2], X[3], X[4]);
}

int main(int argc, char* argv[])
{
    example1();
    example2();
    printf("\n");

    return 0;
}
