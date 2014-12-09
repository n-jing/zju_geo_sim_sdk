#ifndef HJ_MATH_FUNC_OPERATION_H_
#define HJ_MATH_FUNC_OPERATION_H_


#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>

#include <hjlib/sparse/fast_AAT.h>

#include "math_func.h"

namespace hj { namespace math_func {

template <typename FUNC>
class nil_ptr
{
};

template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR>
class math_func_op_const : public math_func_t<VAL_TYPE, INT_TYPE>
{
public:
	math_func_op_const(const math_func &src)
		:src_(src) {
	}
	math_func_op_const(PTR<const math_func> src)
		:src_(*src), own_(src) {
	}
	virtual size_t nx(void) const { return src_.nx(); }
	virtual size_t nf(void) const { return src_.nf(); }

	virtual func_ctx *new_ctx(const void *x) const { return src_.new_ctx(x); }

protected:
	const math_func &src_;
	PTR<const math_func> own_;
};

/*
template <template <typename VAL_TYPE> class OP>
class op_traits
{
public:
	template <typename VAL_TYPE, typename INT_TYPE>
	static int jac(const math_func &func, const VAL_TYPE *x, VAL_TYPE *val, INT_TYPE *ptr, INT_TYPE *idx, func_ctx *ctx, VAL_TYPE scalar) {
		return func.jac(x, val, ptr, idx, ctx);
	}
};

#define HJ_FUNC_OP_ON_JAC(OP)					\
template <>						\
class op_traits<OP>		\
{										\
public:											\
    template <typename VAL_TYPE, typename INT_TYPE>							\
	static int jac(const math_func &func, const VAL_TYPE *x, VAL_TYPE *val, INT_TYPE *ptr, INT_TYPE *idx, func_ctx *ctx, VAL_TYPE scalar) { \
		func.jac(x, val, ptr, idx, ctx);								\
		for(size_t i = 0; i < func.jac_nnz(); ++i)						\
			val[ptr[0]+i] = OP<VAL_TYPE>()(val[ptr[0]+i], scalar); \
		return 0; \
    }	  \
};

HJ_FUNC_OP_ON_JAC(std::multiplies)
HJ_FUNC_OP_ON_JAC(std::divides)

template <template <typename VAL_TYPE> class OP, typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR>
class scalar_op : public math_func_op_const<VAL_TYPE, INT_TYPE, PTR>
{
public:
	typedef VAL_TYPE val_type;
	typedef INT_TYPE int_type;

	scalar_op(const math_func &src, val_type scalar)
		:math_func_op_const<val_type, int_type, PTR>(src), scalar_(scalar) {
	}
	scalar_op(PTR<const math_func> src, val_type scalar)
		:math_func_op_const<val_type, int_type, PTR>(src), scalar_(scalar) {
	}
	virtual int val(const val_type *x, val_type *f, func_ctx *ctx = 0) const {
		math_func_op_const<val_type, int_type, PTR>::src_.val(x, f, ctx);
		for(size_t i = 0; i < math_func_op_const<val_type, int_type, PTR>::nf(); ++i)
			f[i] = OP<val_type>()(f[i], scalar_);
		return 0;
	}
	virtual int jac(const val_type *x, val_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
		return op_traits<OP>::jac(math_func_op_const<val_type, int_type, PTR>::src_, x, val, ptr, idx, ctx, scalar_);
	}
protected:
	val_type scalar_;
};

#define HJ_SCALAR_OP_REF(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE> \
math_func_t<VAL_TYPE, INT_TYPE> *operator OP (const math_func_t<VAL_TYPE, INT_TYPE> &func, VAL_TYPE scalar) { \
	return new scalar_op<NAME, VAL_TYPE, INT_TYPE, nil_ptr>(func, scalar);		\
}
HJ_SCALAR_OP_REF(+, std::plus)
HJ_SCALAR_OP_REF(-, std::minus)
HJ_SCALAR_OP_REF(*, std::multiplies)
HJ_SCALAR_OP_REF(/, std::divides)

#define HJ_SCALAR_OP_OWN(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR> \
math_func_t<VAL_TYPE, INT_TYPE> *operator OP (PTR<const math_func_t<VAL_TYPE, INT_TYPE> > func, VAL_TYPE scalar) { \
	return new scalar_op<NAME, VAL_TYPE, INT_TYPE, PTR>(func, scalar); \
}
HJ_SCALAR_OP_OWN(+, std::plus)
HJ_SCALAR_OP_OWN(-, std::minus)
HJ_SCALAR_OP_OWN(*, std::multiplies)
HJ_SCALAR_OP_OWN(/, std::divides)

template <template <typename VAL_TYPE> class OP, typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR>
class vector_op : public math_func_op_const<VAL_TYPE, INT_TYPE, PTR>
{
public:
	typedef VAL_TYPE val_type;
	typedef INT_TYPE int_type;

	vector_op(const math_func &src, const val_type *vec)
		:math_func_op_const<val_type, int_type, PTR>(src), vec_(vec) {
	}

	virtual int val(const val_type *x, val_type *f, func_ctx *ctx = 0) const {
		math_func_op_const<val_type, int_type, PTR>::src_.val(x, f, ctx);
		for(size_t i = 0; i < math_func_op_const<val_type, int_type, PTR>::nf(); ++i)
			f[i] = OP<val_type>()(f[i], vec_[i]);
		return 0;
	}
	virtual int jac(const val_type *x, val_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
		math_func_op_const<val_type, int_type, PTR>::src_.jac(x, val, ptr, idx, ctx);
		return 0;
	}
protected:
	const val_type *vec_;
};

#define HJ_FUNC_VECTOR_OP_REF(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE> \
math_func_t<VAL_TYPE, INT_TYPE> *operator OP (const math_func_t<VAL_TYPE, INT_TYPE> &func, const VAL_TYPE *vec) { \
	return new vector_op<NAME, VAL_TYPE, INT_TYPE, nil_ptr>(func, vec);			\
}

HJ_FUNC_VECTOR_OP_REF(+, std::plus)
HJ_FUNC_VECTOR_OP_REF(-, std::minus)

#define HJ_FUNC_VECTOR_OP_OWN(OP, NAME) \
template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR> \
math_func_t<VAL_TYPE, INT_TYPE> *operator OP (PTR<const math_func_t<VAL_TYPE, INT_TYPE> > func, const VAL_TYPE *vec) { \
	return new vector_op<NAME, VAL_TYPE, INT_TYPE, PTR>(func, vec); \
}
HJ_FUNC_VECTOR_OP_OWN(+, std::plus)
HJ_FUNC_VECTOR_OP_OWN(-, std::minus)
*/
// [XXX] -> XXX
template <typename VAL_TYPE, typename INT_TYPE, class CON>
class fcat : public math_func_t<VAL_TYPE, INT_TYPE>
{
public:
	typedef VAL_TYPE val_type;
	typedef INT_TYPE int_type;
  typedef CON con_type;
	typedef typename con_type::const_iterator const_iterator;

	fcat(std::shared_ptr<const con_type> funcs)
		:funcs_(funcs) {
    init();
	}

	virtual size_t nx(void) const {
		return (*(funcs_->begin()))->nx();
	}
	virtual size_t nf(void) const {
		return nf_;
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

	virtual func_ctx *new_ctx(const val_type *x) const {
		std::auto_ptr<catenated_func_ctx> ctx(new catenated_func_ctx);
		for(const_iterator i = funcs_->begin();
			i != funcs_->end(); ++i) {
			ctx->push_back((*i)->new_ctx(x));
		}
		return ctx.release();
	}
	virtual int eval(size_t k, const val_type *x, const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const {
		if(ctx) {
      assert(0); // TODO: not tested
			catenated_func_ctx *ctxs = dynamic_cast<catenated_func_ctx *>(ctx);
			assert(ctxs && ctxs->size() == funcs_->size());
			typename catenated_func_ctx::container::iterator itr = ctxs->begin();
			for(const_iterator i = funcs_->begin();
				i != funcs_->end(); ++i, ++itr) {
				if((*i)->eval(k, x, cv, *itr))
					return 1;
			}
		}
		else {
      size_t i;
#if HJ_MATH_FUNC_USE_OMP
#pragma omp parallel for private(i)
#endif
      for(i = 0; i < f_base_.size(); ++i) {
        coo2val_t<val_type, int_type> cvi = cv;
        cvi.f_base() += f_base_[i];
        (*funcs_)[i]->eval(k, x, cvi, 0);
      }
      return 0;
    }
  }
	virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
		if(ctx) {
      assert(0); // TODO: not tested
			// catenated_func_ctx *ctxs = dynamic_cast<catenated_func_ctx *>(ctx);
			// assert(ctxs && ctxs->size() == funcs_->size());
			// typename catenated_func_ctx::container::iterator itr = ctxs->begin();
			// for(const_iterator i = funcs_->begin(); i != funcs_->end(); ++i, ++itr) {
			// 	if((*i)->patt(k, x, coo, *itr))
			// 		return 1;
			// }
		}
		else {
      coo_l2g l2gi = l2g;
      for( const_iterator i = funcs_->begin(); i != funcs_->end(); ++i) {
        assert((*i)->nnz(k) >= 0); // dense will not go here
        (*i)->patt(k, cs, l2gi, 0);
        l2gi.f_base() += (*i)->nf();
      }
      return 0;
		}
	}
	virtual size_t nnz(size_t k) const {
    const_iterator i = funcs_->begin();
    const size_t first = (*i)->nnz(k); // same sparse/density to the first one.
    if(first == -2) return -2;
    size_t nnz = first;
    for(++i; i != funcs_->end(); ++i) {
      const size_t nnzi = (*i)->nnz(k);
      if(nnzi == -2) return -2;
      if(first == -1 && nnzi != -1)
        return -2;
      if(first != -1) {
        if(nnzi == -1)
          return -2;
        else
          nnz += nnzi;
      }
    }
    return nnz;
  }
protected:
	void init(void) {
		assert(funcs_->size());
    f_base_.resize(funcs_->size());
    f_base_[0] = 0;
    for(size_t i = 1; i < f_base_.size(); ++i)
      f_base_[i] = f_base_[i-1]+(*funcs_)[i-1]->nf();
		nf_ = f_base_.back()+funcs_->back()->nf();
	}
	size_t nf_;
  std::vector<size_t> f_base_;
  std::shared_ptr<const con_type> funcs_;
};

template <typename VAL_TYPE, typename INT_TYPE, class CON>
math_func_t<VAL_TYPE, INT_TYPE> *new_fcat(const std::shared_ptr<CON> &con)
{
	return new fcat<VAL_TYPE, INT_TYPE, CON>(std::shared_ptr<const CON>(con));
}

namespace for_hj_sparse {
template <typename VAL_TYPE, typename INT_TYPE>
class ptr_csc
{
public:
  typedef VAL_TYPE val_type;
  typedef INT_TYPE int_type;

  ptr_csc(size_t nr, size_t nc, size_t nnz,
          const int_type *ptr, const int_type *idx, val_type *val)
    :nr_(nr), nc_(nc), nnz_(nnz), ptr_(ptr), idx_(idx), val_(val) {}

  size_t size(int dim) const { return (dim == 1)?nr_:nc_; }
  size_t nnz(void) const { return nnz_; }
  const int_type* ptr(void) const { return ptr_; }
  const int_type* idx(void) const { return idx_; }
  const val_type* val(void) const { return val_; }
  val_type* val(void) { return val_; }
private:
  const int_type *ptr_, *idx_;
  val_type *val_;
  size_t nr_, nc_, nnz_;
};
}

//! NOTICE: w will be squared \|w*f\|^2
// DS ->DDS
template <typename VAL_TYPE, typename INT_TYPE>
class sumsqr : public math_func_t<VAL_TYPE, INT_TYPE>
{
public:
	typedef VAL_TYPE val_type;
	typedef INT_TYPE int_type;

	sumsqr(const std::shared_ptr<const math_func> &f,
         const std::shared_ptr<const std::vector<VAL_TYPE> > w
         = std::shared_ptr<const std::vector<VAL_TYPE> >(0))
		:f_(f), w_(w) {
    init();
	}

	virtual size_t nx(void) const {
		return f_->nx();
	}
	virtual size_t nf(void) const {
		return 1;
	}
  virtual int eval(size_t k, const val_type *x,
                   const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    matrix<val_type> r = zeros<val_type>(f_->nf(), 1);
    if(f_->eval(0, x, coo2val(*cp_[0], &r[0])))
      return __LINE__;
    if(w_.get()) {
      for(size_t i = 0; w_.get() && i < r.size(); ++i)
        r[i] *= (*w_)[i];
    }
    if(k == 0) {
      int_type c[] = {0};
      cv[c] += dot(r, r);
      return 0;
    }

    matrix<val_type> JT_val = zeros<double>(hj::sparse::nnz(JT_), 1);
    f_->eval(1, x, coo2val(*cp_[1], &JT_val[0]));
    if(w_.get()) {
      for(size_t fi = 0; fi < JT_.size(2); ++fi) {
        for(size_t nzi = JT_.ptr()[fi]; nzi < JT_.ptr()[fi+1]; ++nzi)
          JT_val[nzi] *= (*w_)[fi];
      }
    }
    for_hj_sparse::ptr_csc<val_type, int_type> JT
      (JT_.size(1), JT_.size(2), JT.nnz(), &JT_.ptr()[0], &JT_.idx()[0], &JT_val[0]);
    if(k == 1) {
      matrix<val_type> JTr = zeros<double>(JT.size(1), 1);
      hj::sparse::mv(false, JT, r, JTr);
      for(int_type i = 0; i < JTr.size(); ++i) {// not worth omp
        int_type c[] = {0, i};
        cv[c] += JTr[i]*2;
      }
      return 0;
    }
    if(k == 2) {
      matrix<val_type> JTJ_val = zeros<double>(hj::sparse::nnz(H_), 1);
      for_hj_sparse::ptr_csc<val_type, int_type> JTJ
        (H_.size(1), H_.size(2), H_.nnz(), &H_.ptr()[0], &H_.idx()[0], &JTJ_val[0]);
      fast_AAT(JT, JTJ, true);
      int_type ci;
      for(ci = 0; ci < H_.size(2); ++ci) { // not worth omp
        for(size_t nzi = H_.ptr()[ci]; nzi < H_.ptr()[ci+1]; ++nzi) {
          int_type c[3] = {0, ci, H_.idx()[nzi]};
          cv[c] += JTJ_val[nzi]*2;
        }
      }
    }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
    if(k == 2) {
      size_t cooi = 0;
      for(size_t ci = 0; ci < H_.size(2); ++ci) {
        for(size_t nzi = H_.ptr()[ci]; nzi < H_.ptr()[ci+1]; ++nzi, ++cooi) {
          int_type c[3] = {0, int_type(ci), H_.idx()[nzi]};
          l2g.add(cs, c);
        }
      }
    }
    return 0;
  }
  virtual size_t nnz(size_t k) const {
    if(k == 0)
      return -1;
    if(k == 1)
      return -1;
    if(k == 2) // Gauss-Newton
      return hj::sparse::nnz(H_);
  }
protected:
	void init(void) {
    using namespace hj::sparse;
    using namespace zjucad::matrix;
    for(size_t k = 0; k < 2; ++k)
      cp_[k].reset(hj::math_func::patt<int_type>(*f_, k));
    JT_.resize(f_->nx(), f_->nf(), cp_[1]->nnz());
    JT_.val()(colon()) = 1;
    coo2csc(*cp_[1], &JT_.ptr()[0], &JT_.idx()[0]);
    assert(hj::sparse::nnz(H_) == 0);
    AAT<map_by_sorted_vector>(JT_, H_);
  }
  std::shared_ptr<const math_func> f_;
  std::shared_ptr<const std::vector<VAL_TYPE> > w_;
  hj::sparse::csc<val_type, int_type> JT_, H_;
  std::shared_ptr<coo_pat<int_type> > cp_[2];
};

// DDD -> DSS
template <typename VAL_TYPE, typename INT_TYPE>
class xmap : public math_func_t<VAL_TYPE, INT_TYPE>
{
public:
	typedef VAL_TYPE val_type;
	typedef INT_TYPE int_type;

	xmap(const std::shared_ptr<const math_func> &f, const std::vector<INT_TYPE> &idx, size_t nx)
		:f_(f), idx_(idx), nx_(nx) {
    init();
	}

	virtual size_t nx(void) const {
		return nx_;
	}
	virtual size_t nf(void) const {
		return f_->nf();
	}
  virtual int eval(size_t k, const val_type *x,
                   const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    std::vector<val_type> dx(idx_.size());
    for(size_t i = 0; i < dx.size(); ++i)
      dx[i] = x[idx_[i]];
    std::vector<val_type> v(cp_[k]->nnz(), 0);
    if(f_->eval(k, &dx[0], coo2val(*cp_[k], &v[0])))
      return __LINE__;
    if(k == 0) {
      for(int32_t fi = 0, vi = 0; fi < nf(); ++fi) {
        int32_t c[] = {fi};
        cv[c] += v[vi];
      }
    }
    if(k == 1) {
      for(int32_t fi = 0, vi = 0; fi < nf(); ++fi) {
        for(int32_t xi = 0; xi < idx_.size(); ++xi, ++vi) {
          int32_t c[] = {fi, idx_[xi]};
          cv[c] += v[vi];
        }
      }
    }
    if(k == 2) {
      for(int32_t fi = 0, vi = 0; fi < nf(); ++fi) {
        for(int32_t xi0 = 0; xi0 < idx_.size(); ++xi0) {
          for(int32_t xi1 = 0; xi1 < idx_.size(); ++xi1, ++vi) {
            int32_t c[] = {fi, idx_[xi0], idx_[xi1]};
            cv[c] += v[vi];
          }
        }
      }
    }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
    if(k == 1) {
      for(int32_t fi = 0; fi < nf(); ++fi) {
        for(int32_t xi = 0; xi < idx_.size(); ++xi) {
          int32_t c[] = {fi, idx_[xi]};
          l2g.add(cs, c);
        }
      }
    }
    if(k == 2) {
      for(int32_t fi = 0; fi < nf(); ++fi) {
        for(int32_t xi0 = 0; xi0 < idx_.size(); ++xi0) {
          for(int32_t xi1 = 0; xi1 < idx_.size(); ++xi1) {
            int32_t c[] = {fi, idx_[xi0], idx_[xi1]};
            l2g.add(cs, c);
          }
        }
      }
    }
    return 0;
  }
  virtual size_t nnz(size_t k) const {
    assert(f_->nnz(k) == -1);
    if(k == 0)
      return -1;
    if(k == 1)
      return idx_.size();
    if(k == 2)
      return idx_.size()*idx_.size();
    return -2;
  }
protected:
	void init(void) {
    for(size_t k = 0; k < 3; ++k) { // 0, 1, 2 only
      cp_[k].reset(hj::math_func::patt<int_type>(*f_, k));
    }
  }
  std::shared_ptr<const math_func> f_;
  std::shared_ptr<coo_pat<int_type> > cp_[3];
  const std::vector<int_type> idx_;
  size_t nx_;
};

//! NOTICE: w will not be squared sum w_i*f_i
// -> DDS
template <typename VAL_TYPE, typename INT_TYPE>
class sum : public math_func_t<VAL_TYPE, INT_TYPE>
{
public:
	typedef VAL_TYPE val_type;
	typedef INT_TYPE int_type;

	sum(const std::shared_ptr<const math_func> &f,
         const std::shared_ptr<const std::vector<VAL_TYPE> > w
         = std::shared_ptr<const std::vector<VAL_TYPE> >(0))
		:f_(f), w_(w) {
    init();
	}

	virtual size_t nx(void) const {
		return f_->nx();
	}
	virtual size_t nf(void) const {
		return 1;
	}
  virtual int eval(size_t k, const val_type *x,
                   const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    if(k == 0) {
      matrix<val_type> r = zeros<val_type>(cp_[0]->nnz(), 1);
      if(f_->eval(0, x, coo2val(*cp_[0], &r[0])))
        return __LINE__;
      if(w_.get()) {
        for(size_t i = 0; w_.get() && i < r.size(); ++i)
          r[i] *= (*w_)[i];
      }
      int_type c[] = {0};
      cv[c] += zjucad::matrix::sum(r);
      return 0;
    }
    if(k == 1) {
      matrix<val_type> g = zeros<val_type>(cp_[1]->nnz(), 1);
      if(f_->eval(1, x, coo2val(*cp_[1], &g[0])))
        return __LINE__;
      for(size_t i = 0; i < g.size(); ++i) {
        int_type c[] = {0, cache_[1][i*2+1]};
        cv[c] += g[i]*(w_?(*w_)[cache_[1][i*2+0]]:1);
      }
      return 0;
    }
    if(k == 2) {
      matrix<val_type> h = zeros<val_type>(cp_[2]->nnz(), 1);
      if(f_->eval(2, x, coo2val(*cp_[2], &h[0])))
        return __LINE__;
      for(size_t i = 0; i < h.size(); ++i) {
        int_type c[] = {0, cache_[2][i*3+1], cache_[2][i*3+2]};
        cv[c] += h[i]*(w_?(*w_)[cache_[2][i*3+0]]:1);
      }
      return 0;
    }
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
    if(k == 2) {
      int c2[2];
      for(size_t ci = 0; ci < cp2_->nnz(); ++ci) {
        (*cp2_)(ci, c2);
        int_type c[3] = {0, c2[0], c2[1]};
        l2g.add(cs, c);
      }
    }
    return 0;
  }
  virtual size_t nnz(size_t k) const {
    if(k == 0)
      return -1;
    if(k == 1)
      return -1;
    if(k == 2)
      return cp2_->nnz();
  }
protected:
	void init(void) {
    for(size_t k = 0; k < 3; ++k) {
      cp_[k].reset(hj::math_func::patt<int_type>(*f_, k));
    }
    coo_set<int_type> cs(2, nx(), nx(), cp_[2]->nnz());
    coo_l2g l2g;
    int_type c[3];
    for(size_t i = 0; i < cp_[2]->nnz(); ++i) {
      (*cp_[2])(i, c);
      l2g.add(cs, c+1);
    }
    cp2_.reset(cs.get_coo_pat());
    cache_.resize(3);
    for(int k = 0; k < 3; ++k) {
      size_t nnzk = cp_[k]->nnz();
      std::vector<int_type> coo(k+1);
      cache_[k].reserve(nnzk*coo.size());
      for(size_t i = 0; i < nnzk; ++i) {
        (*cp_[k])(i, &coo[0]);
        for(size_t d = 0; d < coo.size(); ++d)
          cache_[k].push_back(coo[d]);
      }
    }
  }
  std::shared_ptr<const math_func> f_;
  std::shared_ptr<const std::vector<VAL_TYPE> > w_;
  std::shared_ptr<coo_pat<int_type> > cp_[3], cp2_;
  std::vector<std::vector<int_type> > cache_;
};

// sub math_func has dense g and hes, typical use:
// sum x_map<square<dense_fun> >
// DDD -> DDD
template <typename VAL_TYPE, typename INT_TYPE>
class square : public math_func_t<VAL_TYPE, INT_TYPE>
{
public:
	typedef VAL_TYPE val_type;
	typedef INT_TYPE int_type;

	square(const std::shared_ptr<const math_func> &f)
		:f_(f) {
    init();
	}

	virtual size_t nx(void) const {
		return f_->nx();
	}
	virtual size_t nf(void) const {
		return f_->nf();
	}
  virtual int eval(size_t k, const val_type *x,
                   const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    matrix<val_type> r = zeros<val_type>(cp_[0]->nnz(), 1);
    if(f_->eval(0, x, coo2val(*cp_[0], &r[0])))
      return __LINE__;
    if(k == 0) {
      for(int_type i = 0; i < f_->nf(); ++i) {
        int_type c[] = {i};
        cv[c] += r[i]*r[i];
      }
      return 0;
    }
    matrix<val_type> g = zeros<val_type>(cp_[1]->nnz(), 1);
    if(f_->eval(1, x, coo2val(*cp_[1], &g[0])))
      return __LINE__;
    if(k == 1) {
      size_t gi = 0;
      for(int32_t fi = 0; fi < nf(); ++fi, gi+=nx()) {
        for(int32_t xi = 0; xi < nx(); ++xi) {
          int32_t c[] = {fi, xi};
          cv[c] += g[gi+xi]*r[fi]*2;
        }
      }
      return 0;
    }
    // no H
    if(k == 2) { // not worth omp
      size_t gi = 0;
      for(int32_t fi = 0; fi < nf(); ++fi, gi += nx()) {
        for(int32_t xi0 = 0; xi0 < nx(); ++xi0) {
          for(int32_t xi1 = 0; xi1 < nx(); ++xi1) {
            int32_t c[] = {fi, xi0, xi1};
            cv[c] += g[gi+xi0]*g[gi+xi1]*2;
          }
        }
      }
      return 0;
    }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
    return 0;
  }
  virtual size_t nnz(size_t k) const {
    assert(f_->nnz(k) == -1);
    return -1;
  }
protected:
	void init(void) {
    for(size_t k = 0; k < 3; ++k) { // 0, 1, 2 only
      cp_[k].reset(hj::math_func::patt<int_type>(*f_, k));
      if(!cp_[k].get())
        throw f_.get();
    }
  }
  std::shared_ptr<const math_func> f_;
  std::shared_ptr<coo_pat<int_type> > cp_[3];
  const std::vector<int_type> idx_;
  size_t nx_;
};

}}

#endif
