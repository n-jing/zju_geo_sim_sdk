#ifndef HJ_FUNCTION_OPERATION_H_
#define HJ_FUNCTION_OPERATION_H_

#include "function.h"

namespace hj { namespace function {

template <typename FUNC>
class nil_ptr
{
};

template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR>
class function_op_const : public function_t<VAL_TYPE, INT_TYPE>
{
public:
	function_op_const(const function &src)
		:src_(src) {
	}
	function_op_const(PTR<const function> src)
		:src_(*src), own_(src) {
	}
	virtual size_t dim_of_x(void) const { return src_.dim_of_x(); }
	virtual size_t dim_of_f(void) const { return src_.dim_of_f(); }

	virtual func_ctx *new_ctx(const void *x) const { return src_.new_ctx(x); }

	virtual size_t jac_nnz(void) const {
		return src_.jac_nnz();
	}
protected:
	const function &src_;
	PTR<const function> own_;
};

template <template <typename VAL_TYPE> class OP>
class op_traits
{
public:
	template <typename VAL_TYPE, typename INT_TYPE>
	static int jac(const function &func, const VAL_TYPE *x, VAL_TYPE *val, INT_TYPE *ptr, INT_TYPE *idx, func_ctx *ctx, VAL_TYPE scalar) {
		return func.jac(x, val, ptr, idx, ctx);
	}
};

#define HJ_FUNC_OP_ON_JAC(OP)					\
template <>						\
class op_traits<OP>		\
{										\
public:											\
    template <typename VAL_TYPE, typename INT_TYPE>							\
	static int jac(const function &func, const VAL_TYPE *x, VAL_TYPE *val, INT_TYPE *ptr, INT_TYPE *idx, func_ctx *ctx, VAL_TYPE scalar) { \
		func.jac(x, val, ptr, idx, ctx);								\
		for(size_t i = 0; i < func.jac_nnz(); ++i)						\
			val[ptr[0]+i] = OP<VAL_TYPE>()(val[ptr[0]+i], scalar); \
		return 0; \
    }	  \
};

HJ_FUNC_OP_ON_JAC(std::multiplies)
HJ_FUNC_OP_ON_JAC(std::divides)

template <template <typename VAL_TYPE> class OP, typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR>
class scalar_op : public function_op_const<VAL_TYPE, INT_TYPE, PTR>
{
public:
	typedef VAL_TYPE value_type;
	typedef INT_TYPE int_type;

	scalar_op(const function &src, value_type scalar)
		:function_op_const<value_type, int_type, PTR>(src), scalar_(scalar) {
	}
	scalar_op(PTR<const function> src, value_type scalar)
		:function_op_const<value_type, int_type, PTR>(src), scalar_(scalar) {
	}
	virtual int val(const value_type *x, value_type *f, func_ctx *ctx = 0) const {
		function_op_const<value_type, int_type, PTR>::src_.val(x, f, ctx);
		for(size_t i = 0; i < function_op_const<value_type, int_type, PTR>::dim_of_f(); ++i)
			f[i] = OP<value_type>()(f[i], scalar_);
		return 0;
	}
	virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
		return op_traits<OP>::jac(function_op_const<value_type, int_type, PTR>::src_, x, val, ptr, idx, ctx, scalar_);
	}
protected:
	value_type scalar_;
};

#define HJ_SCALAR_OP_REF(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE> \
function_t<VAL_TYPE, INT_TYPE> *operator OP (const function_t<VAL_TYPE, INT_TYPE> &func, VAL_TYPE scalar) { \
	return new scalar_op<NAME, VAL_TYPE, INT_TYPE, nil_ptr>(func, scalar);		\
}
HJ_SCALAR_OP_REF(+, std::plus)
HJ_SCALAR_OP_REF(-, std::minus)
HJ_SCALAR_OP_REF(*, std::multiplies)
HJ_SCALAR_OP_REF(/, std::divides)

#define HJ_SCALAR_OP_OWN(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR> \
function_t<VAL_TYPE, INT_TYPE> *operator OP (PTR<const function_t<VAL_TYPE, INT_TYPE> > func, VAL_TYPE scalar) { \
	return new scalar_op<NAME, VAL_TYPE, INT_TYPE, PTR>(func, scalar); \
}
HJ_SCALAR_OP_OWN(+, std::plus)
HJ_SCALAR_OP_OWN(-, std::minus)
HJ_SCALAR_OP_OWN(*, std::multiplies)
HJ_SCALAR_OP_OWN(/, std::divides)

template <template <typename VAL_TYPE> class OP, typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR>
class vector_op : public function_op_const<VAL_TYPE, INT_TYPE, PTR>
{
public:
	typedef VAL_TYPE value_type;
	typedef INT_TYPE int_type;

	vector_op(const function &src, const value_type *vec)
		:function_op_const<value_type, int_type, PTR>(src), vec_(vec) {
	}

	virtual int val(const value_type *x, value_type *f, func_ctx *ctx = 0) const {
		function_op_const<value_type, int_type, PTR>::src_.val(x, f, ctx);
		for(size_t i = 0; i < function_op_const<value_type, int_type, PTR>::dim_of_f(); ++i)
			f[i] = OP<value_type>()(f[i], vec_[i]);
		return 0;
	}
	virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
		function_op_const<value_type, int_type, PTR>::src_.jac(x, val, ptr, idx, ctx);
		return 0;
	}
protected:
	const value_type *vec_;
};

#define HJ_FUNC_VECTOR_OP_REF(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE> \
function_t<VAL_TYPE, INT_TYPE> *operator OP (const function_t<VAL_TYPE, INT_TYPE> &func, const VAL_TYPE *vec) { \
	return new vector_op<NAME, VAL_TYPE, INT_TYPE, nil_ptr>(func, vec);			\
}

HJ_FUNC_VECTOR_OP_REF(+, std::plus)
HJ_FUNC_VECTOR_OP_REF(-, std::minus)

#define HJ_FUNC_VECTOR_OP_OWN(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR> \
function_t<VAL_TYPE, INT_TYPE> *operator OP (PTR<const function_t<VAL_TYPE, INT_TYPE> > func, const VAL_TYPE *vec) { \
	return new vector_op<NAME, VAL_TYPE, INT_TYPE, PTR>(func, vec); \
}
HJ_FUNC_VECTOR_OP_OWN(+, std::plus)
HJ_FUNC_VECTOR_OP_OWN(-, std::minus)

template <typename VAL_TYPE, typename INT_TYPE, template <typename CON> class PTR, typename CON>
class catenated_function : public function_t<VAL_TYPE, INT_TYPE>
{
public:
	typedef VAL_TYPE value_type;
	typedef INT_TYPE int_type;
	typedef typename CON::const_iterator const_iterator;

	catenated_function(const CON &funcs)
		:funcs_(funcs) {
		init();
	}
	catenated_function(PTR<const CON> funcs)
		:funcs_(*funcs), own_(funcs) {
		init();
	}

	virtual size_t dim_of_x(void) const {
		return (*funcs_.begin())->dim_of_x();
	}
	virtual size_t dim_of_f(void) const {
		return dim_of_f_;
	}

	class catenated_func_ctx : public func_ctx
	{
	public:
		typedef std::vector<func_ctx *> container;

		void push_back(func_ctx *ctx) {
			ctxs_.push_back(ctx);
		}
		size_t size(void) const {
			return ctxs_.size();
		}
		container::iterator begin(void) {
			return ctxs_.begin();
		}
		virtual ~catenated_func_ctx() {
			for(container::iterator i = ctxs_.begin();
				i != ctxs_.end(); ++i)
				delete *i;
		}
	protected:
		container ctxs_;
	};

	virtual func_ctx *new_ctx(const value_type *x) const {
		std::auto_ptr<catenated_func_ctx> ctx(new catenated_func_ctx);
		for(const_iterator i = funcs_.begin();
			i != funcs_.end(); ++i) {
			ctx->push_back((*i)->new_ctx(x));
		}
		return ctx.release();
	}
	virtual int val(const value_type *x, value_type *f, func_ctx *ctx = 0) const {
		if(ctx) {
			catenated_func_ctx *ctxs = dynamic_cast<catenated_func_ctx *>(ctx);
			assert(ctxs && ctxs->size() == funcs_.size());
			typename catenated_func_ctx::container::iterator itr = ctxs->begin();
			for(const_iterator i = funcs_.begin();
				i != funcs_.end(); ++i, ++itr) {
				if((*i)->val(x, f, *itr))
					return 1;
				f += (*i)->dim_of_f();
			}
		}
		else {
			for(const_iterator i = funcs_.begin();
				i != funcs_.end(); ++i) {
				if((*i)->val(x, f, 0))
					return 1;
				f += (*i)->dim_of_f();
			}
		}
		return 0;
	}
	virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
		if(ctx) {
			catenated_func_ctx *ctxs = dynamic_cast<catenated_func_ctx *>(ctx);
			assert(ctxs && ctxs->size() == funcs_.size());
			typename catenated_func_ctx::container::iterator itr = ctxs->begin();
			for(const_iterator i = funcs_.begin();
				i != funcs_.end(); ++i, ++itr) {
				if((*i)->jac(x, val, ptr, idx, *itr))
					return 1;
				if(ptr)
					ptr += (*i)->dim_of_f();
			}
		}
		else {
			for(const_iterator i = funcs_.begin();
				i != funcs_.end(); ++i) {
				if((*i)->jac(x, val, ptr, idx, 0))
					return 1;
				if(ptr)
					ptr += (*i)->dim_of_f();
			}
		}
		return 0;
	}
	virtual size_t jac_nnz(void) const {
		return jac_nnz_;
	}
protected:
	void init(void) {
		assert(funcs_.size());
		dim_of_f_ = 0;
		jac_nnz_ = 0;
		for(const_iterator i = funcs_.begin();
			i != funcs_.end(); ++i) {
			dim_of_f_ += (*i)->dim_of_f();
			jac_nnz_ += (*i)->jac_nnz();
		}
	}
	size_t dim_of_f_, jac_nnz_;
	const CON &funcs_;
	PTR<const CON> own_;
};

template <typename VAL_TYPE, typename INT_TYPE, typename CON>
function_t<VAL_TYPE, INT_TYPE> *new_catenated_function(const CON &con)
{
	return new catenated_function<VAL_TYPE, INT_TYPE, nil_ptr, CON>(con);
}

template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR, typename CON>
function_t<VAL_TYPE, INT_TYPE> *new_catenated_function(PTR<CON> &con)
{
	return new catenated_function<VAL_TYPE, INT_TYPE, PTR, CON>(PTR<const CON>(con));
}

}}

#endif
