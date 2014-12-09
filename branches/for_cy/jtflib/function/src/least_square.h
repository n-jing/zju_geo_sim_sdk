#ifndef HJ_LEAST_SQUARE_H_
#define HJ_LEAST_SQUARE_H_

#include <hjlib/function/function.h>

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>

#include "../include/function.h"

namespace jtf {
  namespace function {

//! assume the dim and pattern of f is constant
class least_square_Gauss_Newton : public jtf::function::functionN1_t<double,int32_t>
{
public:
  least_square_Gauss_Newton(const hj::function::function_t<double, int32_t> &f);
  least_square_Gauss_Newton(const std::shared_ptr<const hj::function::function_t<double, int32_t> > own);
  virtual size_t dim(void) const;
  virtual int val(const double *x, double &v);

  //! @param nnz=0 means dense, g=0 or h=0 means query nnz
  virtual int gra(const double *x, double *g);
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx);
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1);
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    return -1;
  }
protected:
  const hj::function::function_t<double, int32_t> &f_;
  const std::shared_ptr<const hj::function::function_t<double, int32_t> > own_;
  mutable zjucad::matrix::matrix<double> r_;
  hj::sparse::csc<double, int32_t> JT_, hes_;
  std::set<size_t> idx_;
  bool is_JT_sorted_;
  bool is_JT_created;
};

  }}
#endif
