#ifndef JTF_MATH_H
#define JTF_MATH_H

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>

namespace jtf{
  namespace math{
    void null_basis(const zjucad::matrix::matrix<double> &A,
                    zjucad::matrix::matrix<double> &z,
                    const double threshold);

//    void null_basis(const hj::sparse::csc<double,int32_t> &A,
//                    hj::sparse::csc<double,int32_t> &z,
//                    const double threshold);
  }
}
#endif
