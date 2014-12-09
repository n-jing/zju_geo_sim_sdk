#include "../include/optimizer.h"
#include "Gauss-Newton.h"
#include "alglib.h"
#include "libLBFGS.h"
#include "HLBFGS.h"

#include <string>
#include <iostream>

using namespace hj::function;
using namespace std;

namespace zjucad {

int optimize(
	const hj::function::function &f,
	zjucad::matrix::matrix<double> &x,
	zjucad::matrix::matrix<double> &residual,
	boost::property_tree::ptree& pt
	)
{
	residual.resize(f.dim_of_f());
	pt.put("iter.desc","number of max iteration");

	pt.put("package.desc","numerical package:<alglib,libLBFGS,hj,HLBFGS>");
	const string package = pt.get<string>("package.value");
	if(package == "alglib") {
	  pt.put("alg.desc","algorithm <non-linear-cg,lbfgs>");
		const string alg = pt.get<string>("alg.value", "algorithm <non-linear-cg,lbfgs> ");
		if(pt.get<string>("alg.value") == "non-linear-cg") {
			int rtn = alglib_min(f, x, residual, pt);
			cerr << "alglib_min returns: " << rtn << endl;
		}
		else if(pt.get<string>("alg.value") == "lbfgs") {
			int rtn = alglib_min(f, x, residual, pt);
			cerr << "alglib_min returns: " << rtn << endl;
		}
	}
	else if(package == "libLBFGS") {
		int rtn = libLBFGS(f, x, residual, pt);
		cerr << "libLBFGS returns: " << rtn << endl;
	}
	else if(package == "hj") {
		const string alg = pt.get<string>("alg.value", "More");
		if(alg == "More") {
			More_Gauss_Newton opt(f);
			opt.iterate(&x[0], f, &residual[0], pt);
		}
	}
	else if(package == "HLBFGS") {
		HLBFGS(f, x, residual, pt);
	}
	else {
		cerr << "# no optimizer found." << endl;
	}
	return 0;
}

}
