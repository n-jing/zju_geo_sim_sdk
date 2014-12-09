#ifndef HJ_FUNCTION_AUX_H_
#define HJ_FUNCTION_AUX_H_

#include "func_abi.h"

namespace hj { namespace function {

extern "C" {
HJ_FUNCTION_API int hj_calc_jac(const function_in_c *func, const function *this_, void *vx, void *vjac, void *veps);
HJ_FUNCTION_API double hj_jac_err(const function_in_c *func, const function *this_, void *vx);
}

//! jac is #f*#x
template <typename VAL_TYPE>
int calc_jac(const function &f, const VAL_TYPE *x, VAL_TYPE *jac,
			 VAL_TYPE eps = std::numeric_limits<VAL_TYPE>::epsilon()) {
	function_in_c cf;
	return hj_calc_jac(&cf,
					   &f,
					   reinterpret_cast<void *>(const_cast<VAL_TYPE *>(x)),
					   reinterpret_cast<void *>(jac),
					   reinterpret_cast<void *>(&eps));
}

template <typename VAL_TYPE>
double jac_err(const function &f, const VAL_TYPE *x) {
	function_in_c cf;
	return hj_jac_err(&cf,
                    &f,
                    reinterpret_cast<void *>(const_cast<VAL_TYPE *>(x)));
}

}}

#endif
