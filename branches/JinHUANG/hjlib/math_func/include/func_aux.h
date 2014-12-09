#ifndef HJ_MATH_FUNC_AUX_H_
#define HJ_MATH_FUNC_AUX_H_

#include "math_func.h"

namespace hj { namespace math_func {

//! X ends with -1.
//! assume k-1 is correct, coo is a length 2 integer (f_(k)_i, xi)
int numerical_check(const math_func &f, size_t k, const void *x,
                    void *err, void *coo, const void *eps,
                    size_t parallel_num, const void *X);

template <typename VAL_TYPE, typename INT_TYPE>
VAL_TYPE numerical_check(const math_func &f, size_t k, const VAL_TYPE *x,
                         INT_TYPE *coo, const VAL_TYPE *eps = 0,
                         size_t parallel_num = 8, const INT_TYPE *X = 0) {
  assert(k > 0);
  VAL_TYPE err = 0;
  if(numerical_check(f, k, reinterpret_cast<const void *>(x),
                     reinterpret_cast<void *>(&err), reinterpret_cast<void *>(coo),
                     reinterpret_cast<const void *>(eps), parallel_num, X))
    return -1;
  return err;
  
}

}}

#endif
