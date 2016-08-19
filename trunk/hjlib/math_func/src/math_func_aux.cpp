#include <algorithm>
#include <iostream>
#include <vector>

#include <zjucad/matrix/include/matrix.h>
#include <zjucad/matrix/include/io.h>

#include "../include/func_aux.h"
#include "../include/coo.h"

namespace hj { namespace math_func {

using namespace std;
using namespace zjucad::matrix;

template <typename VAL_TYPE, typename INT_TYPE>
int numerical_check_tt(const math_func &f, size_t k, const VAL_TYPE *x0,
                       VAL_TYPE *err, INT_TYPE *coo, const VAL_TYPE *eps,
                       size_t parallel_num, const INT_TYPE *X)
{
  assert(x0 && err && coo && X);
  unique_ptr<coo_pat<INT_TYPE> > coo_k[2];
  for(size_t i = 0; i < 2; ++i) {
    const size_t nnz_k = f.nnz(k-1+i);
    coo_set<INT_TYPE> s(k+i, f.nf(), f.nx(), nnz_k);
    f.patt(k-1+i, s, coo_l2g());
    coo_k[i].reset(s.get_coo_pat());
  }
  matrix<VAL_TYPE> ref = zeros<VAL_TYPE>(coo_k[1]->nnz(), 1); {
    f.eval(k, x0, coo2val(*coo_k[1], &ref[0]));
  }
  matrix<VAL_TYPE> ct = zeros<VAL_TYPE>(coo_k[0]->nnz(), 1); {
    f.eval(k-1, x0, coo2val(*coo_k[0], &ct[0]));
  }
  vector<matrix<INT_TYPE> > coo0(ct.size());
  size_t i;
#pragma omp parallel for private(i)
  for(i = 0; i < ct.size(); ++i) {
    coo0[i].resize(k+1);
    (*coo_k[0])(i, &coo0[i][0]);
  }

  size_t xn;
  for(xn = 0; X[xn] != -1; ++xn);
  size_t gi;
  vector<pair<VAL_TYPE, vector<INT_TYPE> > > sub_max(parallel_num);
#pragma omp parallel for private(gi)
  for(gi = 0; gi < parallel_num; ++gi) {
    sub_max[gi].first = 0;
    sub_max[gi].second.resize(k+1);
    vector<VAL_TYPE> xgi(f.nx());
    matrix<VAL_TYPE> ct_eps;
    copy(x0, x0+xgi.size(), xgi.begin());
    for(size_t xi = gi; xi < xn; xi+=parallel_num) {
      const VAL_TYPE save = xgi[X[xi]];
      xgi[X[xi]] += *eps;
      ct_eps = -ct;
      f.eval(k-1, &xgi[0], coo2val(*coo_k[0], &ct_eps[0]));
      ct_eps /= *eps;
      for(size_t i = 0; i < ct_eps.size(); ++i) {
        matrix<INT_TYPE> coo_i = coo0[i];
        coo_i[k] = X[xi];
        const size_t ad = (*coo_k[1])(&coo_i[0]);
        VAL_TYPE e = ct_eps[i];
        if(ad != -1)
          e -= ref[ad];
        if(fabs(e) > fabs(sub_max[gi].first)) {
          sub_max[gi].first = e;
          copy(coo_i.begin(), coo_i.end(), sub_max[gi].second.begin());
        }
      }
      xgi[X[xi]] = save;
    }
  }

  typename vector<pair<VAL_TYPE, vector<INT_TYPE> > >::const_iterator max
    = max_element(sub_max.begin(), sub_max.end());
  *err = max->first;
  copy(max->second.begin(), max->second.end(), coo);
  return 0;
}

template <typename VAL_TYPE>
int numerical_check_t(const math_func &f, size_t k, const VAL_TYPE *x,
                      VAL_TYPE *err, void *coo, const VAL_TYPE *eps,
                      size_t parallel_num, const void *X = 0)
{
#define CALL_I(T) \
  if(f.get_int_byte() == sizeof(T)) {                                   \
    VAL_TYPE default_eps = 1e-8;                                        \
    matrix<T> X0;                                                       \
    if(!X) {                                                            \
      X0 = colon(0, f.nx());                                            \
      X0[f.nx()] = -1;                                                \
    }                                                                   \
    return numerical_check_tt                                           \
      (f, k, x, err, reinterpret_cast<T *>(coo),                        \
       (eps?eps:&default_eps), parallel_num,                            \
       reinterpret_cast<const T *>(X?X:&X0[0]));                        \
  }

  CALL_I(int32_t);
  CALL_I(int64_t);
  return __LINE__;
}

int numerical_check(const math_func &f, size_t k, const void *x,
                    void *err, void *coo,
                    const void *eps, size_t parallel_num, const void *X)
{
#define CALL_V(T) \
  if(f.get_value_byte() == sizeof(T))                                 \
    return numerical_check_t(f, k, reinterpret_cast<const T *>(x),    \
      reinterpret_cast<T *>(err), coo,                                \
      reinterpret_cast<const T *>(eps), parallel_num, X);
  CALL_V(double);
  CALL_V(float);
  return __LINE__;
}

}}
