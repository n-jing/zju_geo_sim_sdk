#ifndef HJ_FUNCTION_H_
#define HJ_FUNCTION_H_

#include <stdint.h>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <memory>

#include "config.h"
#include "type.h"

namespace hj { namespace function {

#ifdef __GNUG__
#  define HJ_FUNC_DEPRECATED __attribute__((deprecated))
#endif

#ifdef _MSC_VER
#  define HJ_FUNC_DEPRECATED __declspec(deprecated)
#endif

template <typename T> inline char type2char(void);

// often used as cache
class HJ_FUNC_DEPRECATED func_ctx
{
public:
	virtual ~func_ctx(){}
};

//! generic function without known value and int type to avoid template
class HJ_FUNC_DEPRECATED function
{
public:
	virtual size_t dim_of_x(void) const = 0;
	virtual size_t dim_of_f(void) const = 0;

	template <typename VAL_TYPE>
	func_ctx *new_ctx(const VAL_TYPE *x) const { 
		assert(get_value_type() == type2char<VAL_TYPE>());
		return new_ctx(reinterpret_cast<const void *>(x));
	}

	template <typename VAL_TYPE>
	int val(const VAL_TYPE *x, VAL_TYPE *f, func_ctx *ctx = 0) const {
		assert(get_value_type() == type2char<VAL_TYPE>());
		return val(reinterpret_cast<const void *>(x), reinterpret_cast<void *>(f), ctx);
	}

	template <typename VAL_TYPE, typename INT_TYPE>
	int jac(const VAL_TYPE *x, VAL_TYPE *val, INT_TYPE *ptr = 0, INT_TYPE *idx = 0, func_ctx *ctx = 0) const {
		assert(get_value_type() == type2char<VAL_TYPE>());
		assert(get_int_type() == type2char<INT_TYPE>());
		return jac(reinterpret_cast<const void *>(x), reinterpret_cast<void *>(val),
				   reinterpret_cast<void *>(ptr), reinterpret_cast<void *>(idx), ctx);
	}

	virtual size_t jac_nnz(void) const {
		return dim_of_x()*dim_of_f();
	}

	virtual char get_value_type(void) const = 0;
	virtual char get_int_type(void) const = 0;

	virtual ~function(){}
//protected:
	virtual func_ctx *new_ctx(const void *x) const { return 0; }

	virtual int val(const void *x, void *f, func_ctx *ctx = 0) const = 0;
	//! @brief a |x|*|f| csc matrix with dim_of_nz_x()*|f| nnz
	virtual int jac(const void *x, void *val, void *ptr = 0, void *idx = 0, func_ctx *ctx = 0) const = 0;
};

//! generic function with known value and int type
template <typename VAL_TYPE, typename INT_TYPE>
class HJ_FUNC_DEPRECATED function_t : public function
{
public:
	typedef VAL_TYPE value_type;
	typedef INT_TYPE int_type;

	virtual func_ctx *new_ctx(const value_type *x) const { return 0; }
	virtual int val(const value_type *x, value_type *f, func_ctx *ctx = 0) const = 0;
	virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const = 0;

	virtual char get_value_type(void) const {
		return type2char<value_type>();
	}
	virtual char get_int_type(void) const {
		return type2char<int_type>();
	}

protected:
	virtual func_ctx *new_ctx(const void *x) const {
		return new_ctx(reinterpret_cast<const value_type *>(x));
	}
	virtual int val(const void *x, void *f, func_ctx *ctx = 0) const {
		return val(reinterpret_cast<const value_type *>(x), reinterpret_cast<value_type *>(f), ctx);
	}
	virtual int jac(const void *x, void *val, void *ptr = 0, void *idx = 0, func_ctx *ctx = 0) const {
		return jac(reinterpret_cast<const value_type *>(x), reinterpret_cast<value_type *>(val),
				   reinterpret_cast<int_type *>(ptr), reinterpret_cast<int_type *>(idx), ctx);
	}
};

}}

#include "operation.h"

#endif
