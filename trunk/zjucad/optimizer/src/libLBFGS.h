#ifndef HJ_LIBLBFGS_H_
#define HJ_LIBLBFGS_H_

#include <hjlib/function/function.h>
#include <zjucad/matrix/matrix.h>
#include <boost/property_tree/ptree.hpp>

//! eps : x, f, g
int libLBFGS(const hj::function::function &fun, zjucad::matrix::matrix<double> &x,
			 zjucad::matrix::matrix<double> &residual,
		boost::property_tree::ptree &opts);

#endif
