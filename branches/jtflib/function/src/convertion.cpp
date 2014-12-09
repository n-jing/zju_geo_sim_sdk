#include "../include/function.h"

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>

#include <hjlib/function/function.h>
#include <stdexcept>

namespace jtf{
  namespace function {

    class hj2jtf_func : public jtf::function::functionN1_t<double,int32_t>
    {
    public:
      hj2jtf_func(const std::shared_ptr<
                  const hj::function::function_t<double,int32_t> > src)
        :src_(src){
        if(!src_ || src_->dim_of_f() != 1)
          throw std::logic_error(" input hj_function should be RN->1");
      }
      virtual ~hj2jtf_func(){}

      virtual size_t dim(void) const{
        return src_->dim_of_x();
      }
      virtual int val(const double *x, double &v) {
        double vv = 0;
        src_->val(x, &vv);
        v += vv;
        return 0;
      }
      virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx)
      {
        if(g == 0 && idx == 0) {
            nnz = src_->jac_nnz();
            return 0;
          }

        zjucad::matrix::matrix<int32_t> ptr = zjucad::matrix::zeros<int32_t>(2,1);
        return src_->jac(x, g, &ptr[0], idx);
      }
      virtual int gra(const double *x, double *g){
        zjucad::matrix::itr_matrix<double*> g_m(dim(),1,g);
        zjucad::matrix::matrix<int32_t> ptr = zjucad::matrix::zeros<int32_t>(2,1);
        size_t nnz = src_->jac_nnz();
        zjucad::matrix::matrix<double> g_v = zjucad::matrix::zeros<double>(nnz,1);
        zjucad::matrix::matrix<int32_t> idx = zjucad::matrix::zeros<int32_t>(nnz,1);

        if(src_->jac(x, &g_v[0], &ptr[0], &idx[0]))
          return __LINE__;

        for(size_t i = 0; i < idx.size(); ++i)
          g[idx[i]] += g_v[i];
        return 0;
      }\
      virtual int hes(const double *x, size_t &nnz, size_t &format,
                      double *h, int32_t *ptr, int32_t *idx, double alpha)
      {
        // do not know how to get hes
        if(h == 0 && ptr == 0 && idx == 0)
          nnz = 0;
        return 0;
      }
      virtual int hes_block(const double *x, double *h, double alpha)
      {
        return __LINE__;
      }
    private:
      const std::shared_ptr<const hj::function::function_t<double,int32_t> > src_;
    };


    class jtf2hj_func : public hj::function::function_t<double,int32_t>
    {
    public:
      jtf2hj_func(const std::shared_ptr<
                  const jtf::function::functionN1_t<double,int32_t> > src)
        :src_(src), nnz_(-1){
        if(!src_)
          throw std::invalid_argument("empty jtf::function pointer.");
      }
      virtual ~jtf2hj_func(){}
      virtual size_t dim_of_x()const{ return src_->dim();}
      virtual size_t dim_of_f()const{return 1;}
      virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0)const
      {
        zjucad::matrix::itr_matrix<double*> f0(dim_of_f(),1,f);
        f0 *= 0;
        return const_cast<jtf::function::functionN1_t<double,int32_t>*>(src_.get())->val(x, *f);
      }
      virtual int jac(const double *x, double *val, int32_t *ptr, int32_t *idx,
                      hj::function::func_ctx *ctx = 0) const
      {
        ptr[1] += ptr[0] + jac_nnz();
        size_t nnz = 0;
        zjucad::matrix::itr_matrix<double*> val_m(jac_nnz(),1, val);
        val_m *= 0;
        std::vector<value_type> val_(jac_nnz());
        std::vector<int_type> idx_(jac_nnz());

        const_cast<jtf::function::functionN1_t<double,int32_t>*>(src_.get())->gra(x, nnz, &val_[0], &idx_[0]);

        for(size_t i = 0; i < jac_nnz(); ++i){
            idx[ptr[0]+i] = idx_[i];
            val[ptr[0]+i] = val_[i];
          }
        return 0;
      }
      virtual size_t jac_nnz()const{
        if(nnz_ == -1){
            static zjucad::matrix::matrix<double> fake_x(dim_of_x(),1);
            const_cast<jtf::function::functionN1_t<double,int32_t>*>(src_.get())->gra(&fake_x[0], nnz_, 0,0);
          }
        return nnz_;
      }
    private:
      const std::shared_ptr<const jtf::function::functionN1_t<double,int32_t> > src_;
      mutable size_t nnz_ ;
    };


    jtf::function::functionN1_t<double,int32_t> *
    hj_func_to_jtf_func(const hj::function::function_t<double,int32_t> * func)
    {
      std::unique_ptr<jtf::function::functionN1_t<double,int32_t> > p(
            new hj2jtf_func(std::shared_ptr<const hj::function::function_t<double,int32_t> >(func)));
      return p.release();
    }

    hj::function::function_t<double,int32_t> *
    jtf_func_to_hj_func(const jtf::function::functionN1_t<double,int32_t> *func)
    {
      std::unique_ptr<hj::function::function_t<double,int32_t> > p(
            new jtf2hj_func(std::shared_ptr<const jtf::function::functionN1_t<double,int32_t> >(func)));
      return p.release();
    }
  }
}
