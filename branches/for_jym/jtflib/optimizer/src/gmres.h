#ifndef JTF_GMRES_H
#define JTF_GMRES_H

#include <string>
#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>
#include "precondition.h"

namespace jtf{

  // dense version
  class gmres_solver{
  public:

    gmres_solver(const zjucad::matrix::matrix<double> & A,
                 const zjucad::matrix::matrix<double> & b)
      :A_(A), b_(b), eps_(1e-6){}

    ///
    /// @brief solve Ax=b, dense storge and with out preconditions
    ///        QR factorization is using Given rotations
    /// @param x
    /// @param b
    /// @param m
    /// @param max_iter
    ///
    int solve(double *x, size_t m, size_t max_iter) const;



    int solve2(double *x, size_t m, size_t max_iter,
               const std::shared_ptr<jacobi_precondition> pre)const;

  protected:
  private:
    const zjucad::matrix::matrix<double> &A_;
    const zjucad::matrix::matrix<double> &b_;
    const double eps_;
  };

  class gmres_solver_sparse{
  public:

    gmres_solver_sparse(const hj::sparse::csc<double,int32_t> & A,
                        const zjucad::matrix::matrix<double> & b)
      :A_(A), b_(b), eps_(1e-6){}

    ///
    /// @brief solve Ax=b, dense storge and with out preconditions
    ///        QR factorization is using Given rotations
    /// @param x
    /// @param b
    /// @param m
    /// @param max_iter
    ///
    int solve(double *x, size_t m, size_t max_iter) const;


    int solve2(double *x, size_t m, size_t max_iter,
               const std::shared_ptr<jacobi_precondition_sparse> pre)const;

  protected:
  private:
    const hj::sparse::csc<double,int32_t> & A_;
    const zjucad::matrix::matrix<double> &b_;
    const double eps_;
  };


}
#endif
