#include "precondition.h"
#include <zjucad/matrix/itr_matrix.h>

using namespace zjucad::matrix;

namespace jtf{
  void jacobi_precondition::init()
  {
    M_ = ones<double>(A_.size(1),1);
    const size_t dim = std::min(A_.size(1),A_.size(2));
    for(size_t i = 0; i < dim; ++i)
      M_[i] = weight_/A_(i,i);
  }


  void jacobi_precondition::get_r(const double *x, double *r) const
  {
    itr_matrix<const double*> x_m(A_.size(2),1,x);
    itr_matrix<double*> r_m(b_.size(),1,r);
    r_m = b_-A_*x_m;

    get_b(r);
  }

  void jacobi_precondition::get_b(double *b) const
  {
    assert(b);
    for(size_t i = 0; i < M_.size(); ++i)
      b[i] *= M_[i];
  }

  void jacobi_precondition::get_Ax(const double *x, double *Ax)const
  {
    // M*A*x
    assert(x && Ax);
    itr_matrix<const double*> x_m(A_.size(2),1,x);
    itr_matrix<double*> Ax_m(A_.size(1),1,Ax);
    //A*x
    Ax_m = A_ * x_m;
    get_b(Ax);
  }


  void jacobi_precondition_sparse::init()
  {
    M_ = ones<double>(A_.size(1),1);
    const size_t dim = std::min(A_.size(1),A_.size(2));
    for(size_t i = 0; i < dim; ++i)
      M_[i] = weight_/A_(i,i);
  }


  void jacobi_precondition_sparse::get_r(const double *x, double *r) const
  {
    assert(x && r);
    itr_matrix<double*> r_m(A_.size(1),1,r);
    itr_matrix<const double*> x_m(A_.size(2),1,x);
    r_m = -1*b_;
    hj::sparse::mv(false, A_, x_m, r_m);
    r_m *= -1;
    get_b(r);
  }

  void jacobi_precondition_sparse::get_b(double *b) const
  {
    assert(b);
    for(size_t i = 0; i < M_.size(); ++i)
      b[i] *= M_[i];
  }

  void jacobi_precondition_sparse::get_Ax(const double *x, double *Ax)const
  {
    // M*A*x
    assert(x && Ax);
    itr_matrix<const double*> x_m(A_.size(2),1,x);
    itr_matrix<double*> Ax_m(A_.size(1),1,Ax);
    //A*x
    Ax_m *= 0;
    hj::sparse::mv(false, A_, x_m, Ax);

    get_b(Ax);
  }

}
