#ifndef HJ_LINEAR_SOLVER_H_
#define HJ_LINEAR_SOLVER_H_

#include <boost/property_tree/ptree.hpp>
#include <hjlib/sparse/sparse.h>

class linear_solver
{
public:
  /** 
   * 
   * 
   * @param double 
   * @param A 
   * @param opts "linear_solver/type.value" "<PETsc, direct, iterative>"
   *             "linear_solver/name.value" "<cholmod, umfpack>"
   *             "PETsc/pc.value", "<sor>" for petsc only
   * @return linear_solver pointer
   */
	static linear_solver *create(
		const hj::sparse::csc<double, int32_t> &A,
		boost::property_tree::ptree &opts);

	virtual int solve(const double *b, double *x, size_t rhs, boost::property_tree::ptree &opts) = 0;

	virtual ~linear_solver(){}
};

#endif
