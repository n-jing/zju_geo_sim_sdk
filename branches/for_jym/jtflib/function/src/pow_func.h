//#ifndef POW_FUNC_H
//#define POW_FUNC_H

//#include "../include/function.h"
//#include "../include/func_aux.h"
//#include <stdexcept>
//#include <zjucad/matrix/matrix.h>
//#include <zjucad/matrix/itr_matrix.h>
//#include <map>

//namespace jtf{
//  namespace function{

//    //! @brief: for inequality, -log(f)
//    class pow_func : public jtf::function::functionN1_t<double,int32_t>
//    {
//    public:
//      pow_func(std::shared_ptr<const jtf::function::functionN1_t<val_type,int_type> > own,
//               val_type p)
//        :own_(own), f_(*own), p_(p),format_(-1), nnz_(-1){}
//      pow_func(const jtf::function::functionN1_t<val_type,int_type> &f, val_type p)
//        :f_(f), p_(p),format_(-1), nnz_(-1){}
//      virtual ~pow_func(){}
//      virtual size_t dim(void) const { return f_.dim();}
//      virtual int val(const val_type *x, val_type &v) const
//      {
//        double vv = 0;
//        f_.val(x, vv);
//        v += pow(vv,p_);
//        return 0;
//      }
//      virtual int gra(const val_type *x, val_type *g){
//        double vv = 0;
//        f_.val(x, vv);

//        static zjucad::matrix::matrix<val_type> gra_(dim(),1);
//        if(gra_.size() != dim()) gra_.resize(dim(),1);
//        gra_ *= 0;

//        f_.gra(x, &gra_[0]);

//        zjucad::matrix::itr_matrix<val_type*> g_m(dim(),1, g);
//        g_m += p_ *pow(vv,p_-1) * gra_;
//        return 0;
//      }
//      virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx){
//        if(g == 0 && idx == 0){
//            return f_.gra(x, nnz, g, idx);
//          }
//        if(g != 0 && idx != 0){
//            double vv = 0;
//            f_.val(x, vv);
//            static zjucad::matrix::matrix<val_type> gra_(nnz,1);
//            if(gra_.size() != nnz) gra_.resize(nnz,1);
//            gra_ *= 0;

//            f_.gra(x, nnz, &gra_[0], idx);
//            zjucad::matrix::itr_matrix<val_type*> g_m(nnz,1, g);
//            g_m += p_ * pow(vv, p_-1) * gra_;
//            return 0;
//          }
//      }
//      virtual int hes(const val_type *x, size_t &nnz, size_t &format, val_type *h,
//                      int_type *ptr, int_type *idx, val_type alpha = 1){
//        if(h == 0){
//            if(col_row_.empty()){
//                col_row_.resize(dim());
//                size_t nnz_, format_;
//                f_.hes(x, nnz_, format_, 0,0, 0);

//                std::vector<int_type> idx_(nnz_);
//                std::vector<int_type> ptr_(dim()+1);
//                ptr_[0] = 0;

//                f_.hes(x, nnz, format, 0, &ptr_[0], &idx_[0]);
//                if(format == 1){ //csc
//                    for(size_t i = 1; i < ptr_.size(); ++i){
//                        if(ptr_[i-1] == ptr[i]) continue;
//                        for(size_t t = ptr[i-1]; t < ptr[i]; ++t)
//                          col_row_[i].insert(idx_[t]);
//                      }
//                  }else if(format == 2){ // i j k
//                    for(size_t i = 0; i < nnz; ++i){
//                        col_row_[ptr_[i]].insert(idx_[i]);
//                      }
//                  }

//                size_t jac_nnz;
//                f_.gra(x, jac_nnz,0,0);
//                std::vector<val_type> g(jac_nnz);
//                std::vector<int_type> gra_idx(jac_nnz);
//                f_.gra(x, jac_nnz, &g[0], &gra_idx[0]);
//                for(size_t ji = 0; ji < jac_nnz; ++ji)
//                  for(size_t jj = 0; jj < jac_nnz; ++jj)
//                    col_row_[gra_idx[ji]].insert(gra_idx[jj]);
//              }

//            if(ptr == 0 && idx == 0){
//                nnz = 0;
//                for(const auto & one_col: col_row_){
//                    nnz += one_col.size();
//                  }
//                format = 1;
//                return 0;
//              }else{

//                for(size_t i = 0; i < col_row_.size(); ++i){
//                    ptr[i+1] = ptr[i] + col_row_[i].size();
//                    size_t idx_i = ptr[i];
//                    for(const auto & row_idx : col_row_[i]){
//                        idx[idx_i++] = row_idx;
//                      }
//                  }
//                return 0;
//              }
//          }else{
//            assert(h != 0 && ptr != 0 && idx != 0);

//            size_t nnz_, format_;
//            f_.hes(x, nnz_, format_, 0,0,0);

//            double vv = 0;
//            f_.val(x, vv);

//            if(nnz_ != 0) {

//                static zjucad::matrix::matrix<val_type> h_m_0(nnz_, 1);
//                static zjucad::matrix::matrix<int_type> ptr_m(dim()+1,1),idx_m(nnz_,1);
//                if(h_m_0.size() != nnz_) h_m_0.resize(nnz_,1);
//                if(idx_m.size() != nnz_) idx_m.resize(nnz_,1);
//                h_m_0 *= 0;

//                f_.hes(x, nnz_, format_, 0,&ptr_m[0], &idx_m[0]);
//                f_.hes(x, nnz_, format_, &h_m_0[0], &ptr_m[0], &idx_m[0], alpha);

//                h_m_0 *= p_ * pow(vv,p_-1);

//                if(format_ == 1){
//                    for(int_type i = 0; i < dim(); ++i){
//                        for(int_type r = ptr_m[i]; r < ptr_m[i+1]; ++r){
//                            jtf::function::add_to_csc(h, ptr, idx, i, idx_m[r], h_m_0[r]);
//                          }
//                      }
//                  }else{
//                    for(size_t i = 0; i < nnz_; ++i){
//                        jtf::function::add_to_csc(h, ptr, idx, ptr_m[i], idx_m[i], h_m_0[i]);
//                      }
//                  }
//              }

//            size_t jac_nnz;
//            f_.gra(x, jac_nnz, 0, 0);
//            static zjucad::matrix::matrix<val_type> gra(jac_nnz,1);
//            gra *= 0;
//            static zjucad::matrix::matrix<int_type> gra_idx(jac_nnz,1);
//            gra_idx *= 0;

//            f_.gra(x, jac_nnz, &gra[0], &gra_idx[0]);
//            const val_type a = p_*(p_-1) * pow(vv, p_-2);

//            for(size_t i = 0; i < gra_idx.size(); ++i){
//                for(size_t j = 0; j < gra_idx.size(); ++j){
//                    jtf::function::add_to_csc(h, ptr, idx, gra_idx[i], gra_idx[j],
//                                              gra[i] * gra[j]*a);
//                  }
//              }
//          }
//        return 0;
//      }
//      virtual int hes_block(const val_type *x, val_type *h, val_type alpha = 1) {return -1;}
//      virtual bool is_valid(const val_type *x) const {
//        throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");}
//    protected:
//      double p_;
//      const jtf::function::functionN1_t<val_type,int_type> &f_;
//      const std::shared_ptr<const jtf::function::functionN1_t<val_type,int_type> > own_;
//      size_t format_;
//      size_t nnz_;
//      std::vector<std::set<size_t> > col_row_;
//    };

//  }
//}

//#endif
