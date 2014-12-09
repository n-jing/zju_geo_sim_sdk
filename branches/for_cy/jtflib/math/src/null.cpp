#include "../include/null.h"
#include <stdexcept>
#include <vector>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>


using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;

namespace jtf{
  namespace math{

    inline void spones(zjucad::matrix::matrix<double> &A)
    {
      for(size_t i = 0; i < A.size(); ++i){
          if(fabs(A[i]) > 1e-6)
            A[i] = 1;
        }
    }

    inline void spones(hj::sparse::csc<double,int32_t> &A)
    {
      for(size_t i = 0; i < A.val().size(); ++i){
          if(fabs(A.val()[i]) > 1e-6)
            A.val()[i] = 1;
        }
    }

    inline void sum_matlab(const zjucad::matrix::matrix<double> &A,
                           zjucad::matrix::matrix<double> &r,
                           int type = 1)
    {
      if(type == 1){ // sum each colon
          r.resize(1, A.size(2));
          for(size_t ci = 0; ci < A.size(2); ++ci){
              r[ci] = std::accumulate(A(colon(),ci).begin(), A(colon(),ci).end(),0.0);
            }
        }else if(type == 2){ // sum each row
          r.resize(A.size(1), 1);
          for(size_t ri = 0; ri < A.size(1); ++ri){
              r[ri] = std::accumulate(A(ri,colon()).begin(), A(ri,colon()).end(), 0.0);
            }
        }else{
          throw std::invalid_argument("only support row/col sum");
        }
    }

    inline void sum_matlab(const hj::sparse::csc<double,int32_t> &A,
                           zjucad::matrix::matrix<double> &r,
                           int type = 1)
    {
      if(type == 1){ // sum each colon
          r.resize(1, A.size(2));
          for(size_t ci = 0; ci < A.size(2); ++ci){
              r[ci] = std::accumulate(A.val().begin() + A.ptr()[ci],
                                      A.val().begin() + A.ptr()[ci+1],0.0);
            }
        }else if(type == 2){ // sum each row
          r.resize(A.size(1), 1);
          hj::sparse::csc<double,int32_t> AT;
          hj::sparse::trans(A,AT);
          for(size_t ci = 0; ci < AT.size(2); ++ci){
              r[ci] = std::accumulate(A.val().begin() + A.ptr()[ci],
                                      A.val().begin() + A.ptr()[ci+1],0.0);
            }
        }else{
          throw std::invalid_argument("only support row/col sum");
        }
    }

    inline void sort_matlab(const zjucad::matrix::matrix<double> &r,
                            const int order,// 0: ascend, 1: descend
                            zjucad::matrix::matrix<size_t> & index,
                            zjucad::matrix::matrix<double> * r_ordered = 0)
    {
      vector<pair<double,size_t> > rp(r.size());
      for(size_t i = 0; i < r.size(); ++i){
          rp[i].first = r[i];
          rp[i].second = i;
        }
      if(order == 0)
        sort(rp.begin(), rp.end());
      else
        sort(rp.begin(), rp.end(), std::greater<pair<double,size_t>>());
      index.resize(r.size(),1);
      if(r_ordered) r_ordered->resize(r.size(),1);
      for(size_t i = 0; i < r.size(); ++i){
          index[i] = rp[i].second;
          if(r_ordered){
              (*r_ordered)[i] = rp[i].first;
            }
        }
    }

    inline int nnz_matlab(const zjucad::matrix::matrix<double> &A){
      int nnz = 0;
      for(size_t i = 0; i < A.size(); ++i){
          if(fabs(A[i]) > 1e-6) ++nnz;
        }
      return nnz;
    }

    inline zjucad::matrix::matrix<size_t> find_matlab(const zjucad::matrix::matrix<double> &A)
    {
      matrix<size_t> index(nnz_matlab(A),1);
      for(size_t i = 0,count = 0; i < A.size(); ++i){
          if(fabs(A[i])>1e-6) index[count++] = i;
        }
      return index;
    }

    inline void find_matlab(const zjucad::matrix::matrix<double> &A,
                            matrix<size_t> & index)
    {
      index.resize(nnz_matlab(A),1);
      for(size_t i = 0,count = 0; i < A.size(); ++i){
          if(fabs(A[i])>1e-6) index[count++] = i;
        }
    }

    inline zjucad::matrix::matrix<size_t>
    larger_than(const zjucad::matrix::matrix<double> &A, const double thresold)
    {
      vector<size_t> idx;
      for(size_t i = 0 ; i < A.size(); ++i)
        if(A[i] > thresold) idx.push_back(i);
      matrix<size_t> idx_mat(idx.size(), 1);
      std::copy(idx.begin(), idx.end(), idx_mat.begin());
      return idx_mat;
    }

    template <typename T1>
    inline zjucad::matrix::matrix<size_t>
    no_less_than(const zjucad::matrix::matrix_expression<T1> &A,
                 const double thresold)
    {
      vector<size_t> idx;
      for(size_t i = 0 ; i < A().size(); ++i)
        if(A()[i] > thresold || fabs(A()[i]-thresold) < 1e-6) idx.push_back(i);
      matrix<size_t> idx_mat(idx.size(), 1);
      std::copy(idx.begin(), idx.end(), idx_mat.begin());
      return idx_mat;
    }

    template <typename T1>
    void no_less_than(const zjucad::matrix::matrix_expression<T1> &A,
                      const double thresold,
                      zjucad::matrix::matrix<size_t> &idx)
    {
      vector<size_t> idx_vec;
      for(size_t i = 0 ; i < A().size(); ++i)
        if(A()[i] > thresold || fabs(A()[i]-thresold) < 1e-6) idx_vec.push_back(i);
      idx.resize(idx_vec.size(),1);
      std::copy(idx_vec.begin(), idx_vec.end(), idx.begin());
    }

    pair<double,size_t> min_matlab(const zjucad::matrix::matrix<double> &A)
    {
      pair<double,size_t> min_pair;
      min_pair.second = min_element(A.begin(), A.end()) - A.begin();
      min_pair.first = A[min_pair.second];
      return min_pair;
    }

    zjucad::matrix::matrix<double>& graft_horizontal(
        const zjucad::matrix::matrix<double> &A,
        const zjucad::matrix::matrix<double> &B,
        zjucad::matrix::matrix<double> &C)
    {
      assert(A.size(1) == B.size(1));
      C.resize(A.size(1), A.size(2) + B.size(2));
      C(colon(),colon(0,A.size(2)-1)) = A;
      C(colon(),colon(A.size(2),B.size(2)+A.size(2)-1)) = B;
      return C;
    }

    zjucad::matrix::matrix<double>& graft_vertical(
        const zjucad::matrix::matrix<double> &A,
        const zjucad::matrix::matrix<double> &B,
        zjucad::matrix::matrix<double> &C)
    {
      assert(A.size(2) == B.size(2));
      C.resize(A.size(1)+B.size(1), B.size(2));
      C(colon(0,A.size(1)-1),colon()) = A;
      C(colon(A.size(1), A.size(1)+B.size(1)-1),colon()) = B;
      return C;
    }

    ///////////////////////////////////////////////////////////////////////////

    void null_basis(const zjucad::matrix::matrix<double> &A,
                    zjucad::matrix::matrix<double> &Hz,
                    const double threshold)
    {
      matrix<double> r, A_bkp = A;
      matrix<size_t> t, jnz,jnz_no_less, j;
      spones(A_bkp);
      sum_matlab(A_bkp, r, 2);
      sort_matlab(r,0,t);

      vector<double> H_vec(A.size(2)*A.size(2));
      itr_matrix<double*> H_mat0(A.size(2),A.size(2),&H_vec[0]);

      H_mat0 = eye<double>(A.size(2));
      matrix<double> AT = trans(A);

      matrix<double> BH,BS;
      matrix<double> Hj;
      vector<double> s(A.size(2));

      for(size_t i = 0; i < A.size(1); ++i){
          itr_matrix<double*> H_mat(A.size(2),A.size(2)-i,&H_vec[0]);
          itr_matrix<double*> s_mat(1, H_mat.size(1)-i, &s[0]);

          s_mat = trans(AT(colon(),t[i])) * H_mat;

          if(nnz_matlab(s_mat) == 0) continue;

          find_matlab(s_mat,jnz);

          no_less_than(
                zjucad::matrix::abs(s_mat(jnz)),
                threshold *zjucad::matrix::max(zjucad::matrix::abs(s_mat)), jnz_no_less);
          j = jnz(jnz_no_less);

          Hj = H_mat(colon(),j);
          spones(Hj);
          sum_matlab(Hj,r);
          size_t jj = min_matlab(r).second;

          const size_t jjj = j[jj];

          for(size_t ri = 0; ri < H_mat.size(1); ++ri){
              H_mat(ri,colon()) -= H_mat(ri,jjj) * s_mat /s_mat[jjj];
            }
          for(size_t jump_i = jjj; jump_i < H_mat.size(2)-1; ++jump_i){
              H_mat(colon(), jump_i) = H_mat(colon(),jump_i+1);
              s_mat(0,jump_i) = s_mat(0,jump_i+1);
            }
          H_vec.resize(H_mat.size(1)*(H_mat.size(2)-1));
          s.resize(s.size()-1);
        }
      itr_matrix<double*> H_mat(A.size(2),H_vec.size()/A.size(2),&H_vec[0]);
      Hz  = H_mat;
    }

    //    void null_basis(const hj::sparse::csc<double,int32_t> &A,
    //                    hj::sparse::csc<double,int32_t> &z,
    //                    const double threshold)
    //    {
    //      matrix<double> r;
    //      hj::sparse::csc<double,int32_t> A_bkp = A;
    //      matrix<size_t> t, jnz, j;
    //      spones(A_bkp);
    //      sum_matlab(A_bkp, r, 2);
    //      sort_matlab(r,0,t);

    //      matrix<double> H = eye<double>(A.size(2));

    //      hj::sparse::csc<double,int32_t> AT;
    //      hj::sparse::trans(A,AT);

    //      matrix<double> BH,BS;
    //      matrix<double> Hj,s;
    //      for(size_t i = 0; i < A.size(1); ++i){
    //          s = trans(AT(colon(),t[i])) * H;

    //          if(nnz_matlab(s) == 0) continue;
    //          jnz = find_matlab(s);

    //          j = jnz(no_less_than(
    //                    zjucad::matrix::abs(s(jnz)),
    //                    threshold *zjucad::matrix::max(zjucad::matrix::abs(s))));

    //          Hj = H(colon(),j);

    //          sum_matlab(spones(Hj),r);
    //          size_t jj = min_matlab(r).second;

    //          const size_t jjj = j[jj];

    //          graft_horizontal(H(colon(),colon(0,jjj-1)), H(colon(),colon(jjj+1, H.size(2)-1)),BH);
    //          graft_horizontal(s(0,colon(0,jjj-1)), s(0,colon(jjj+1, s.size(2)-1)),BS);

    //          H = temp(BH - H(colon(), jjj) * BS/s[jjj]);
    //        }

    //    }





  }
}
