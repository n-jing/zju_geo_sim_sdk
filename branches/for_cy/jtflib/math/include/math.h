#ifndef JTFLIB_MATH_H
#define JTFLIB_MATH_H

#include <cmath>
#include <limits>
#include <zjucad/matrix/matrix.h>

namespace jtf{
  namespace math{
#define numerical_precision 1e-8

    constexpr double My_PI()
    {
      return 3.1415926535897932384626;
    }

    template <typename T >
    T safe_acos(T d){
      return std::acos(std::min(std::max(d, -1.0),1.0));
    }

    template <typename T >
    T safe_asin(T d){
      return std::asin(std::min(std::max(d, -1.0),1.0));
    }

    template <typename T>
    T safe_cot(T d){
      return std::tan(My_PI()*0.5-d);
    }

    template<typename T>
    T cos2deg(const T & cos_)
    {
      return safe_acos(cos_) * 180.0 / My_PI();
    }

    template <typename T>
    bool is_nan(T d){
      return d != d;
    }

    template <typename T>
    bool is_inf(T d){
      return std::fabs(d) == std::numeric_limits<T>::infinity();
    }

    template <typename T>
    void erase_nan_inf(T &d){
      if(is_nan(d)) d = 0;
      if(is_inf(d)) d = 0;
    }

    template <typename T>
    T get_sign(T d){
      if(std::fabs(d) < numerical_precision) return 1.0;
      return (d > 0?1.0:-1.0);
    }

    template <typename T>
    void invert_2dmatrix(zjucad::matrix::matrix_expression<T> &A,
                       double EPS =1e-18)
    {
        using namespace zjucad::matrix;
        assert(A().size(1) == A().size(2));
        assert(A().size(1) == 2);
        assert(A().size(1) > 0);
        if(A().size(1) == 1){
            A()(0,0) = 1.0/A()(0,0);
        }else if(A().size(1) == 2){
            const double det = A()(0,0)*A()(1,1) - A()(0,1)*A()(1,0);
            std::swap(A()(0,0),A()(1,1));
            A()(0,1) *= -1;
            A()(1,0) *= -1;
            A() /= det;
        }
    }

  }
}
#endif
