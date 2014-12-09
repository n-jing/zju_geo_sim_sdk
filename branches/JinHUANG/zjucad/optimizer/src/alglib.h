#ifndef HJ_ALGLIB_H_
#define HJ_ALGLIB_H_

#include <hjlib/function/function.h>
#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>

//! eps : x, f, g
int alglib_min(const hj::function::function &fun, zjucad::matrix::matrix<double> &x,
			   zjucad::matrix::matrix<double> &residual,
		//hj::arg_opts &opts);
		boost::property_tree::ptree &pt);

#endif
