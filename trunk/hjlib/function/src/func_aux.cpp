#include "../include/func_aux.h"

#include <memory>

#include <zjucad/matrix/include/matrix.h>
#include <zjucad/matrix/include/itr_matrix.h>

#include <hjlib/sparse/include/sparse.h>

#include <iostream>

using namespace zjucad::matrix;
using namespace hj::function;
using namespace std;

namespace hj { namespace function {

template <typename VAL_TYPE>
int calc_jac_imp(const function_abi &f, const VAL_TYPE *x, VAL_TYPE *jac, VAL_TYPE eps)
{
	const size_t sz[2] = {f.dim_of_x(), f.dim_of_f()};
	VAL_TYPE *px = const_cast<VAL_TYPE *>(x);
	matrix<VAL_TYPE> v0(sz[1]), v1(sz[1]);
	itr_matrix<VAL_TYPE *> mjac(sz[1], sz[0], jac);
	auto_ptr<func_ctx_abi> ct_ctx(f.new_ctx(x));
	f.val((const void *)(x), (void *)(&v0[0]), ct_ctx.get());
	VAL_TYPE tmp;
	eps = sqrt(eps);
	for(size_t xi = 0; xi < sz[0]; ++xi) {
		tmp = x[xi];
		px[xi] += eps;
		auto_ptr<func_ctx_abi> eps_ctx(f.new_ctx(x));
		f.val(x, &v1[0], eps_ctx.get());
		mjac(colon(), xi) = (v1-v0)/eps;
		px[xi] = tmp;
	}
	return 0;
}

#define CALC_JAC(VAL_TYPE)											\
		return calc_jac_imp(f,										\
							reinterpret_cast<const VAL_TYPE *>(vx),	\
							reinterpret_cast<VAL_TYPE *>(vjac),		\
							*reinterpret_cast<VAL_TYPE *>(veps));

#define IF_CALL_JAC(VAL_TYPE)					\
	else if(f.get_value_type() == type2char<VAL_TYPE>()) CALC_JAC(VAL_TYPE)

int hj_calc_jac(const function_in_c *func, const function *this_, void *vx, void *vjac, void *veps)
{
	function_abi f(*func, *this_);
	assert(vx && vjac && veps);
	if(0) {}
	IF_CALL_JAC(double)
	IF_CALL_JAC(float)
	return 1;
}

template <typename VAL_TYPE, typename INT_TYPE>
double jac_err_imp(const function_abi &f, const VAL_TYPE *x)
{
  matrix<VAL_TYPE> J(f.dim_of_f(), f.dim_of_x());
  calc_jac_imp<VAL_TYPE>(f, x, &J[0], std::numeric_limits<VAL_TYPE>::epsilon());
  hj::sparse::csc<VAL_TYPE, INT_TYPE> JT(f.dim_of_x(), f.dim_of_f(), f.jac_nnz());
  f.jac((const void *)(x), (void*)(&JT.val()[0]), (void *)(&JT.ptr()[0]), (void *)(&JT.idx()[0]));
  for(size_t ci = 0; ci < JT.size(2); ++ci) {
    for(size_t nzi = JT.ptr()[ci]; nzi < JT.ptr()[ci+1]; ++nzi) {
      J(ci, JT.idx()[nzi]) -= JT.val()[nzi];
    }
  }
  return max(fabs(J));
}

#define ELSE_IF(VAL_TYPE, INT_TYPE, func)  \
    else if(f.get_value_type() == type2char<VAL_TYPE>() \
            && f.get_int_type() == type2char<INT_TYPE>()) {\
      return func<VAL_TYPE, INT_TYPE>(f, reinterpret_cast<const VAL_TYPE *>(vx)); \
    }

double hj_jac_err(const function_in_c *func, const function *this_, void *vx)
{
	function_abi f(*func, *this_);
	assert(vx);
  if(0) {}
  ELSE_IF(double, int32_t, jac_err_imp)
  ELSE_IF(double, int64_t, jac_err_imp)
  ELSE_IF(float, int32_t, jac_err_imp)
  ELSE_IF(float, int64_t, jac_err_imp)

  return 1;
}

}}
