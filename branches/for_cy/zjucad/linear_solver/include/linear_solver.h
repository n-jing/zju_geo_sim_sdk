#ifndef HJ_LINEAR_SOLVER_H_
#define HJ_LINEAR_SOLVER_H_

#include <boost/property_tree/ptree.hpp>

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
					const double* val, const int32_t* idx,
					const int32_t* ptr, const size_t nnz,
					const size_t row,	const size_t col,
					boost::property_tree::ptree & opts);

	virtual int solve(const double *b, double *x, size_t rhs, boost::property_tree::ptree &opts) = 0;

	virtual ~linear_solver(){}
};

#endif
