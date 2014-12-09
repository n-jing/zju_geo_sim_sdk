#ifndef HJ_MATH_FUNC_H_
#define HJ_MATH_FUNC_H_

#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <memory>

#include "config.h"
#include "coo.h"

namespace hj { namespace math_func {

// often used as cache
class func_ctx
{
public:
	virtual ~func_ctx(){}
};

//! generic function without known value and int type to avoid template
class math_func
{
public:
  virtual size_t nx(void) const = 0;
  virtual size_t nf(void) const = 0;

	template <typename VAL_TYPE>
	func_ctx *new_ctx(const VAL_TYPE *x) const { 
		assert(get_value_byte() == sizeof(VAL_TYPE));
		return new_ctx(reinterpret_cast<const void *>(x));
	}

  template <typename VAL_TYPE, typename INT_TYPE>
	int eval(size_t k, const VAL_TYPE *x, const coo2val_t<VAL_TYPE, INT_TYPE> &cv,
           func_ctx *ctx = 0) const {
		assert(get_value_byte() == sizeof(VAL_TYPE));
		return eval(k, reinterpret_cast<const void *>(x),
                reinterpret_cast<const void *>(&cv), ctx);
	}
  template <typename INT_TYPE>
	int patt(size_t k, coo_set<INT_TYPE> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
		assert(get_int_byte() == sizeof(INT_TYPE));
		return patt(k, reinterpret_cast<void *>(&cs), reinterpret_cast<const void *>(&l2g), ctx);
  }

  //! @param rtn: -1 for dense, -2 for err, other value for nnz
  virtual size_t nnz(size_t k) const = 0;

	virtual char get_value_byte(void) const = 0;
	virtual char get_int_byte(void) const = 0;

  virtual int eval(size_t k, const void *x, const void *cv, func_ctx *ctx = 0) const = 0;
  virtual int patt(size_t k, void *cs, const void *l2g, func_ctx *ctx = 0) const = 0;

  virtual ~math_func(){}
};

template <typename INT>
coo_pat<INT> *patt(const math_func &f, size_t k)
{
  using namespace std;
  size_t nnz = f.nnz(k);
  if(nnz == -2) return 0;
  coo_set<INT> s(k+1, f.nf(), f.nx(), nnz);
  int rtn = f.patt(k, s, coo_l2g());
  if(rtn) return 0;
  return s.get_coo_pat();
}

//! generic function with known value and int type
template <typename VAL_TYPE, typename INT_TYPE>
class math_func_t : public math_func
{
public:
	typedef VAL_TYPE value_type;
	typedef INT_TYPE int_type;

	virtual func_ctx *new_ctx(const value_type *x) const { return 0; }

	virtual int eval(size_t k, const value_type *x, const coo2val_t<value_type, int_type> &cv,
                   func_ctx *ctx = 0) const = 0;
	virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const = 0;

	virtual char get_value_byte(void) const {
		return sizeof(value_type);
	}
	virtual char get_int_byte(void) const {
		return sizeof(int_type);
	}

protected:
	virtual func_ctx *new_ctx(const void *x) const {
		return new_ctx(reinterpret_cast<const value_type *>(x));
	}
	virtual int eval(size_t k, const void *x, const void *cv, func_ctx *ctx = 0) const {
		return eval(k, reinterpret_cast<const value_type *>(x),
                *reinterpret_cast<const coo2val_t<value_type, int_type> *>(cv), ctx);
	}
	virtual int patt(size_t k, void *cs, const void *l2g, func_ctx *ctx = 0) const {
		return patt(k, *reinterpret_cast<coo_set<int_type> *>(cs),
      *reinterpret_cast<const coo_l2g *>(l2g), ctx);
	}
};

}}

#endif
