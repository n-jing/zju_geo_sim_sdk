#ifndef HJ_FUNCTION_ABI_H_
#define HJ_FUNCTION_ABI_H_

#include "config.h"
#include "function.h"

#include <limits>

namespace hj { namespace function {

inline char function_get_value_type(const function &this_) {
	return this_.get_value_type();
}

inline char function_get_int_type(const function &this_) {
	return this_.get_int_type();
}

inline size_t function_dim_of_x(const function &this_) {
	return this_.dim_of_x();
}

inline size_t function_dim_of_f(const function &this_) {
	return this_.dim_of_f();
}
		
inline func_ctx *function_new_ctx(const function &this_, const void *x) {
	return this_.new_ctx(x);
}

inline void function_delete_ctx(func_ctx *ctx) {
	delete ctx;
}

inline int function_val(const function &this_, const void *x, void *f, func_ctx *ctx) {
	return this_.val(x, f, ctx);
}

inline int function_jac(const function &this_, const void *x, void *val,
						void *ptr, void *idx, func_ctx *ctx) {
	return this_.jac(x, val, ptr, idx, ctx);
}

inline size_t function_jac_nnz(const function &this_) {
	return this_.jac_nnz();
}

struct function_in_c
{
	function_in_c() {
		get_value_type = hj::function::function_get_value_type;
		get_int_type = hj::function::function_get_int_type;
		dim_of_x = hj::function::function_dim_of_x;
		dim_of_f = hj::function::function_dim_of_f;
		new_ctx = hj::function::function_new_ctx;
		delete_ctx = hj::function::function_delete_ctx;
		val = hj::function::function_val;
		jac = hj::function::function_jac;
		jac_nnz = hj::function::function_jac_nnz;
	}
	char (*get_int_type)(const function &this_);
	char (*get_value_type)(const function &this_);
	size_t (*dim_of_x)(const function &this_);
	size_t (*dim_of_f)(const function &this_);
	func_ctx *(*new_ctx)(const function &this_, const void *x);
	void (*delete_ctx)(func_ctx *ctx);
	int (*val)(const function &this_, const void *x, void *f, func_ctx *ctx);
	int (*jac)(const function &this_, const void *x, void *val,
			   void *ptr, void *idx, func_ctx *ctx);
	size_t (*jac_nnz)(const function &this_);
};

class func_ctx_abi : public func_ctx
{
public:
	func_ctx_abi(const function_in_c &abi, func_ctx *ctx)
		:abi_(abi), ctx_(ctx) {
	}
	func_ctx *get(void) {
		return ctx_;
	}
	~func_ctx_abi() {
		abi_.delete_ctx(ctx_);
	}
protected:
	const function_in_c &abi_;
	func_ctx *ctx_;
};

class function_abi : public function
{
public:
	function_abi(const function_in_c &abi, const function &func_)
		:abi_(abi), this_(func_) {
	}
	char get_value_type(void) const {
		return abi_.get_value_type(this_);
	}
	char get_int_type(void) const {
		return abi_.get_int_type(this_);
	}
	size_t dim_of_x(void) const {
		return abi_.dim_of_x(this_);
	}
	size_t dim_of_f(void) const {
		return abi_.dim_of_f(this_);
	}
	func_ctx_abi *new_ctx(const void *x) const {
		return new func_ctx_abi(abi_, abi_.new_ctx(this_, x));
	}
	int val(const void *x, void *f, func_ctx *ctx = 0) const {
    if(ctx)
      return abi_.val(this_, x, f, dynamic_cast<func_ctx_abi *>(ctx)->get());
    return abi_.val(this_, x, f, 0);
	}
	int jac(const void *x, void *val,
			void *ptr, void *idx, func_ctx *ctx = 0) const {
    if(ctx)
      return abi_.jac(this_, x, val, ptr, idx, dynamic_cast<func_ctx_abi *>(ctx)->get());
    return abi_.jac(this_, x, val, ptr, idx, 0);
	}
	size_t jac_nnz(void) const {
		return abi_.jac_nnz(this_);
	}
protected:
	const function_in_c &abi_;
	const function &this_;
};

}}

#endif
