#ifndef ZJUCAD__OPTIMIZER_H_
#define ZJUCAD_OPTIMIZER_H_

#include <boost/property_tree/ptree.hpp>
#include <hjlib/function/function.h>
#include <zjucad/matrix/matrix.h>

namespace zjucad {

int optimize(
	const hj::function::function &f,
	zjucad::matrix::matrix<double> &x,
	zjucad::matrix::matrix<double> &residual,
	boost::property_tree::ptree &pt
	);

}

#endif
