#include "../include/function.h"

//double validate_jac(function &f, const double *x);

// double validate_jac(function &f, const double *x)
// {
// 	const itr_matrix<const double *> x0(f.dim_of_x(), 1, x);
// 	matrix<double> x1, f0(f.dim_of_f(), 1), f1 = f0;
// 	f.x_changed();
// 	f.val(&x0[0], &f0[0]);
	
// 	csc<double, function::idx_type> JT(f.dim_of_x(), f.dim_of_f(), f.jac_nnz()), J;
// 	f.jac(&x0[0], &JT.val_[0], &JT.ptr_[0], &JT.idx_[0]);
// 	trans(JT, J);

// 	const double eps = 1e-10;
// 	matrix<double> Ji;
// 	double rtn = 0;
// 	for(size_t xi = 0; xi < f.dim_of_x(); ++xi) {
// 		x1 = x0;
// 		x1[xi] += eps;
// 		f.x_changed();
// 		f.val(&x1[0], &f1[0]);
// 		Ji = (f1-f0)/eps;
// 		for(ptrdiff_t nzi = J.ptr_[xi]; nzi < J.ptr_[xi+1]; ++nzi)
// 			Ji[J.idx_[nzi]] -= J.val_[nzi];
// 		rtn += dot(Ji, Ji);
// 	}
// 	return rtn;
// }

