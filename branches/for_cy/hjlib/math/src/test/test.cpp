#include <iostream>
using namespace std;

#include <zjucad/matrix/matrix.h>	// must before f2c.h in vc

/*
extern "C" {
#include <CLAPACK/F2CLIBS/f2c.h>
#include <ATLAS/cblas.h>
#include <ATLAS/clapack.h>
#include <CLAPACK/clapack.h>
#undef min
#undef max
}
*/
#include "blas_lapack.h"
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/io.h>
using namespace zjucad::matrix;

#include "polar.h"
using namespace hj;

int main(int argc, char *argv[])
{
	polar3d polar;
	double err;
	for(int i = 0; i < 10000; ++i) {
		matrix<double> A = rand<double>(3), R(3, 3), tmp;
		if(polar(A, R, 2) < 0)
			cout << "polar error." << endl;
		tmp = trans(R)*R-eye<double>(3);
		err = norm(tmp);
		tmp = R*trans(R)-eye<double>(3);
		err+= norm(tmp);
		if(err > 2e-7) {
			cout << "orthogonal error" << err << endl << R*trans(R) << trans(R)*R;
			return 1;
		}
		matrix<double> tmpAR[2] = {A, R};
		double detAR[2];
		if((detAR[1] = det(tmpAR[1])) < 0) {
			if((detAR[0] = det(tmpAR[0])) > 0) {
				cout << "det A R not same." << endl;
				return 3;
			}
			cout << "det error" << A << R << endl << detAR[0] << endl << detAR[1] << endl;
			return 2;
		}
	}
	return 0;
}

