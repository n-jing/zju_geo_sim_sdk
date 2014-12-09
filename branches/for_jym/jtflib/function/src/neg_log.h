#ifndef NEG_LOG_H
#define NEG_LOG_H

#include "../include/function.h"
#include "../include/func_aux.h"
#include <hjlib/sparse/fast_AAT.h>
#include <stdexcept>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/sparse/sparse.h>

namespace jtf{
  namespace function{

    //! @brief: for inequality, -log(f)
    class neg_log_func : public jtf::function::functionN1_t<double,int32_t>
    {
    public:
      neg_log_func(std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > own)
        :own_(own), f_(*own_), format_(-1), nnz_(-1){}
      neg_log_func(jtf::function::functionN1_t<double,int32_t> &f)
        :f_(f), format_(-1), nnz_(-1){}
      virtual ~neg_log_func(){}
      virtual size_t dim(void) const { return f_.dim();}
      virtual int val(const double *x, double &v)
      {
        double vv = 0;
        f_.val(x, vv);
        v += -log(vv);
        return 0;
      }
      virtual int gra(const double *x, double *g){
        double vv = 0;
        f_.val(x, vv);

        f_.gra(x, g);
        zjucad::matrix::itr_matrix<double*> g_m(dim(),1, g);
        g_m /= -vv;
        return 0;
      }
      virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx){
        if(g == 0 && idx == 0){
            return f_.gra(x, nnz, g, idx);
          }
        if(g != 0 && idx != 0){
            double vv = 0;
            f_.val(x, vv);
            f_.gra(x, nnz, g, idx);
            zjucad::matrix::itr_matrix<double*> g_m(nnz,1, g);
            g_m /= -vv;
            return 0;
          }
      }
      virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                      int32_t *ptr, int32_t *idx, double alpha = 1){
        if(h == 0){
            if(col_row_.empty()){
                col_row_.resize(dim());
                size_t nnz_, format_;
                f_.hes(x, nnz_, format_, 0,0, 0);

                std::vector<int32_t> idx_(nnz_);
                std::vector<int32_t> ptr_(dim()+1);
                ptr_[0] = 0;

                f_.hes(x, nnz, format, 0, &ptr_[0], &idx_[0]);
                if(format == 1){ //csc
                    for(size_t i = 1; i < ptr_.size(); ++i){
                        if(ptr_[i-1] == ptr[i]) continue;
                        for(size_t t = ptr[i-1]; t < ptr[i]; ++t)
                          col_row_[i].insert(idx_[t]);
                      }
                  }else if(format == 2){ // i j k
                    for(size_t i = 0; i < nnz; ++i){
                        col_row_[ptr_[i]].insert(idx_[i]);
                      }
                  }

                size_t jac_nnz;
                f_.gra(x, jac_nnz,0,0);
                std::vector<double> g(jac_nnz);
                std::vector<int32_t> gra_idx(jac_nnz);
                f_.gra(x, jac_nnz, &g[0], &gra_idx[0]);
                for(size_t ji = 0; ji < jac_nnz; ++ji)
                  for(size_t jj = 0; jj < jac_nnz; ++jj)
                    col_row_[gra_idx[ji]].insert(gra_idx[jj]);
              }

            if(ptr == 0 && idx == 0){
                nnz = 0;
                for(const auto & one_col: col_row_){
                    nnz += one_col.size();
                  }
                format = 1;
                return 0;
              }else{
                for(size_t i = 0; i < col_row_.size(); ++i){
                    ptr[i+1] = ptr[i] + col_row_[i].size();
                    //for(size_t ii = 0; ii < col_row_[i].size(); ++ii){
                    size_t it = 0;
                    for(const auto & row_idx : col_row_[i]){
                        idx[ptr[i] + it++] = row_idx;
                      }
                  }
                return 0;
              }
          }else{
            assert(h != 0 && ptr != 0 && idx != 0);
            f_.hes(x, nnz, format, h, ptr, idx, alpha);
            double vv = 0;
            f_.val(x, vv);
            zjucad::matrix::itr_matrix<double*> h_m(nnz, 1, h);
            h_m /= -vv;
            size_t jac_nnz;
            f_.gra(x, jac_nnz, 0, 0);
            zjucad::matrix::matrix<double> gra = zjucad::matrix::zeros<double>(jac_nnz,1);
            zjucad::matrix::matrix<int32_t> gra_idx = zjucad::matrix::zeros<int32_t>(jac_nnz,1);

            f_.gra(x, jac_nnz, &gra[0], &gra_idx[0]);
            for(size_t i = 0; i < gra_idx.size(); ++i){
                for(size_t j = 0; j < gra_idx.size(); ++j)
                  jtf::function::add_to_csc(h, ptr, idx, gra_idx[i],gra_idx[j],
                                            gra[i] * gra[j]/(vv*vv));
              }
          }
        return 0;
      }
      virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
      virtual bool is_valid(const double *x) const {
        throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");}
    protected:
      //bool is_interior_point(double v) const;
      jtf::function::functionN1_t<double,int32_t> &f_;
      std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > own_;
      size_t format_;
      size_t nnz_;
      std::vector<std::set<size_t> > col_row_;
    };


    class neg_log_from_hj_func : public jtf::function::functionN1_t<double,int32_t>
    {
    public:
      neg_log_from_hj_func(
          std::shared_ptr<const hj::function::function_t<double,int32_t> > own,
          const zjucad::matrix::matrix<double> &weight)
        :own_(own), f_(*own_), format_(-1), nnz_(-1),weight_(weight){
        src_val_.resize(own_->dim_of_f());
        src_jac_.resize(own_->dim_of_x(), own_->dim_of_f(), own_->jac_nnz());

      }

      virtual ~neg_log_from_hj_func(){}
      virtual size_t dim(void) const { return own_->dim_of_x();}
      virtual int val(const double *x, double &v)
      {
        own_->val(x, &src_val_[0]);
        for(size_t fi = 0; fi < src_val_.size(); ++fi) {
            if(!is_interior_point(src_val_[fi])) continue;
            //      src_val_[fi] = NEG_LOG_MIN;
            double t = -log(src_val_[fi])*weight_[fi];
            //erase_nan(t);
            v += t;
          }
        return 0;
      }
      virtual int gra(const double *x, double *g){
        own_->val(x, &src_val_[0]);
        own_->jac(x, &src_jac_.val()[0], &src_jac_.ptr()[0], &src_jac_.idx()[0]);
        for(size_t fi = 0; fi < src_jac_.size(2); ++fi) {
            if(!is_interior_point(src_val_[fi])) continue;
            //      src_val_[fi] = NEG_LOG_MIN;
            for(size_t nzi = src_jac_.ptr()[fi]; nzi < src_jac_.ptr()[fi+1]; ++nzi) {
                double acc = -src_jac_.val()[nzi]/src_val_[fi]*weight_[fi];
                //  erase_nan(acc);
                g[src_jac_.idx()[nzi]] += acc;
              }
          }
        return 0;
      }
      virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx){
        if(g == 0 && idx == 0){
            own_->jac(x, &src_jac_.val()[0], &src_jac_.ptr()[0], &src_jac_.idx()[0]);
            std::set<size_t> nnz_idx(src_jac_.idx().begin(), src_jac_.idx().end());
            nnz = nnz_idx.size();
            return 0;
          }
        if(g != 0 && idx != 0){
            own_->val(x, &src_val_[0]);
            own_->jac(x, &src_jac_.val()[0], &src_jac_.ptr()[0], &src_jac_.idx()[0]);
            std::map<size_t, double> idx_2_val;
            for(size_t fi = 0; fi < src_jac_.size(2); ++fi) {
                if(!is_interior_point(src_val_[fi])) continue;
                //      src_val_[fi] = NEG_LOG_MIN;
                for(size_t nzi = src_jac_.ptr()[fi]; nzi < src_jac_.ptr()[fi+1]; ++nzi) {
                    double acc = -src_jac_.val()[nzi]/src_val_[fi]*weight_[fi];
                    //  erase_nan(acc);
                    auto it = idx_2_val.find(src_jac_.idx()[nzi]);
                    if(it == idx_2_val.end())
                      idx_2_val[src_jac_.idx()[nzi]] = acc;
                    else
                      it->second += acc;
                  }
              }
            return 0;
          }
      }
      virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                      int32_t *ptr, int32_t *idx, double alpha = 1){
        using namespace zjucad::matrix;
        using namespace std;

        format = 1;
        own_->jac(x, &src_jac_.val()[0], &src_jac_.ptr()[0], &src_jac_.idx()[0]);
        if(hj::sparse::nnz(src_hes_))
          fast_AAT(src_jac_, src_hes_, is_JT_sorted_);
        else {
            is_JT_sorted_ = is_sorted_csc(src_jac_);
            hj::sparse::AAT<hj::sparse::map_by_sorted_vector>(src_jac_, src_hes_);
          }
        if(h == 0 && ptr == 0 && idx == 0) { // query nnz
            nnz = hj::sparse::nnz(src_hes_);
            return 0;
          }
        if(h == 0 && ptr !=0 && idx != 0) { // query patten
            if(nnz < hj::sparse::nnz(src_hes_)) {
                cerr << "nnz err: " << nnz << " " << hj::sparse::nnz(src_hes_) << endl;
                return __LINE__;
              }
            copy(src_hes_.ptr().begin(), src_hes_.ptr().end(), ptr);
            copy(src_hes_.idx().begin(), src_hes_.idx().end(), idx);
            return 0;
          }
        if(h != 0 && ptr !=0 && idx != 0) { // accumulate
            own_->val(x, &src_val_[0]);
            matrix<double> src_g, src_H;
            matrix<int32_t> src_g_idx;
            for(size_t fi = 0; fi < src_jac_.size(2); ++fi) {
                if(!is_interior_point(src_val_[fi])) continue;
                //        src_val_[fi] = NEG_LOG_MIN;
                src_g_idx = src_jac_.idx()(colon(src_jac_.ptr()[fi], src_jac_.ptr()[fi+1]-1));
                src_g = src_jac_.val()(colon(src_jac_.ptr()[fi], src_jac_.ptr()[fi+1]-1));
                const double scale = alpha*weight_[fi]/(src_val_[fi]*src_val_[fi]);
                src_H = src_g*trans(src_g)*scale;
                for(int ci = 0; ci < src_g.size(); ++ci) {
                    for(int ri = 0; ri < src_g.size(); ++ri) {
                        double t = src_H(ri, ci);
                        // erase_nan(t);
                        add_to_csc(h, ptr, idx, src_g_idx[ri], src_g_idx[ci], t);
                      }
                  }
              }
            return 0;
          }
        return __LINE__;
      }
      virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
      virtual bool is_valid(const double *x) const {
        double *srv_val_non_cost = const_cast<double *>(&src_val_[0]);
        own_->val(x, srv_val_non_cost);
        for(size_t fi = 0; fi < src_val_.size(); ++fi) {
            //    if(fabs(weight_[fi]) < 1e-5) continue;
            if(!is_interior_point(src_val_[fi])) return false;
          }
        return true;
      }

      bool is_interior_point(double v) const
      {
        //  const double eps = std::numeric_limits<double>::min();
        const double eps = 1e-8;
        if(v < eps) {
            //    cerr << "# not interior of neg_log_sum_func: " << v << endl;
            return false;
          }
        return true;
      }

    protected:
      const hj::function::function_t<double,int32_t> &f_;
      std::shared_ptr<const hj::function::function_t<double,int32_t> > own_;
      size_t format_;
      size_t nnz_;
      std::vector<std::set<size_t> > col_row_;
      const zjucad::matrix::matrix<double> weight_;
      mutable std::vector<double> src_val_;
      hj::sparse::csc<double, int32_t> src_jac_, src_hes_;
      bool is_JT_sorted_;
    };


  }
}

#endif
