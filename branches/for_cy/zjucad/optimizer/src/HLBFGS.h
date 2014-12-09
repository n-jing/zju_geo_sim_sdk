#ifndef HJ_HBFGS_H_
#define HJ_HBFGS_H_

#include <hjlib/function/function.h>
#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>

int HLBFGS(const hj::function::function &fun, zjucad::matrix::matrix<double> &x,
		   zjucad::matrix::matrix<double> &residual,
		boost::property_tree::ptree &opts);

#endif
