/******************************************************************************
 *
 * File:           qp.h
 *
 * Created:        07/03/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Header file for the Quadratic Programming problem solver
 *
 * Description:    None
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined(_QP_H)
#define _QP_H

/** Solver for a Quadratic Programming problem
 * @param nec Number of equality contraints
 * @param nic Number of inequality constraints
 * @param nwc Number of "weak" constraints
 * @param nv Number of variables
 * @param A Data matrix [nec+nic][nv] for equality and inequality constraints
 * @param B Constant vector [nec+nic] for equality and inequality constraints
 * @param Aw Data matrix [nwc][nv] for weak constraints
 * @param Bw Constant vector [nwc] for weak constraints
 * @param W Vector [nwc] of weights for weak constraints in the objective 
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
int qp_solve(int nec, int nic, int nwc, int nv, double** A, double* B, double** Aw, double* Bw, double* W, double* Xmin, double* Xmax, double* X);

#endif
