extern "C" {
#include <qp.h>
}
#pragma comment (lib, "qp.lib")

#include "qp_wrapper.h"

using namespace zjucad::matrix;

#include <vector>
using namespace std;

int qp_solve(int nec, const matrix<double> &A, const matrix<double> &B,	// (A^T*X+B)_i: = 0 (when i < nec); > 0 (when i >= nec)
			 const matrix<double> &Aw, const matrix<double> &Bw, const matrix<double> &W,	// min\|diag(W)*(Aw^T*X+Bw)\|^2
			 const matrix<double> &Xmin, const matrix<double> &Xmax,	// bound
			 matrix<double> &X)
{
	const int nic = A.size(2)-nec, nwc = Aw.size(2), nv = A.size(1);
	vector<const double*> pA(nec+nic), pAw(nwc);
	for(int i = 0; i < nec+nic; ++i) pA[i] = &A(0, i);
	for(int i = 0; i < nwc; ++i) pAw[i] = &Aw(0, i);
	return qp_solve(nec, nic, nwc, nv,
		const_cast<double **>(&pA[0]), const_cast<double *>(&B[0]),
		const_cast<double **>(&pAw[0]), const_cast<double *>(&Bw[0]), const_cast<double *>(&W[0]),
		const_cast<double *>(&Xmin[0]), const_cast<double *>(&Xmax[0]),
		&X[0]);
}
