#ifndef MATH_FUNC2HJ_FUNC_H
#define MATH_FUNC2HJ_FUNC_H

#include <hjlib/function/function.h>
#include <hjlib/math_func/math_func.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>

enum SMART_PTR{STD_SHARED, RAW};

template <typename T, SMART_PTR SP>
struct math_func_ptr{};

template <typename T>
struct math_func_ptr<T, STD_SHARED>{
  typedef std::shared_ptr<T> math_func_ptr_type;
};

template <typename T>
struct math_func_ptr<T, RAW>{
  typedef T* math_func_ptr_type;
};

template <typename VAL_TYPE, typename INT_TYPE, SMART_PTR SP = STD_SHARED>
class math2hj_func: public hj::function::function_t<VAL_TYPE, INT_TYPE>
{
public:
  typedef hj::math_func::math_func math_func_type;
  typedef typename math_func_ptr<const math_func_type, SP>::math_func_ptr_type mfp_type;

  math2hj_func(mfp_type mptr):math_func_ptr_(mptr){
    for(size_t k = 0; k < 2; ++k){
        cp_[k].reset(hj::math_func::patt<INT_TYPE>(*math_func_ptr_, k));
      }
  }

  virtual size_t dim_of_x(void) const{
    return math_func_ptr_->nx();
  }
  virtual size_t dim_of_f(void) const{
    return math_func_ptr_->nf();
  }

  int val(const VAL_TYPE *x, VAL_TYPE *f, hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::itr_matrix<VAL_TYPE*> f0(dim_of_f(),1,f);
    f0 *= 0;
    math_func_ptr_->eval(0,x, hj::math_func::coo2val(*cp_[0], f));
    return 0;
  }

  int jac(const VAL_TYPE *x, VAL_TYPE *val, INT_TYPE *ptr = 0, INT_TYPE *idx = 0,
          hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::itr_matrix<VAL_TYPE*> val0(jac_nnz(),1,val);
    val0 *= 0;
    hj::math_func::coo2csc(*cp_[1], ptr, idx);
    math_func_ptr_->eval(1, x, hj::math_func::coo2val(*cp_[1], val));

    zjucad::matrix::itr_matrix<const INT_TYPE*> ptr0(dim_of_x()+1,1,ptr);
    zjucad::matrix::itr_matrix<const INT_TYPE*> idx0(jac_nnz(),1,idx);

    return 0;
  }

  virtual size_t jac_nnz(void) const {
    size_t nnz = math_func_ptr_->nnz(1);
    if(nnz == -1) return math_func_ptr_->nx();
    return nnz;
  }

private:
  math2hj_func(){}
  math2hj_func & operator = (const math2hj_func &){}

private:
  mfp_type math_func_ptr_;
  std::shared_ptr<hj::math_func::coo_pat<int32_t> > cp_[3];
};

#endif // MATH_FUNC2HJ_FUNC_H
