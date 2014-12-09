#include <memory>

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>

#include <hjlib/sparse/sparse.h>

#include <iostream>

#include "../include/func_aux.h"
#include "least_square.h"

using namespace zjucad::matrix;
using namespace jtf::function;
using namespace std;

namespace jtf { namespace function {

    jtf::function::functionN1_t<double,int32_t> *
    least_square_warpper(const hj::function::function_t<double,int32_t> & f)
    {
      std::unique_ptr<jtf::function::functionN1_t<double,int32_t> > lf(
            new jtf::function::least_square_Gauss_Newton(f));
      return lf.release();
    }

    jtf::function::functionN1_t<double,int32_t> *
    least_square_warpper(const std::shared_ptr<const hj::function::function_t<double,int32_t> > fptr)
    {
      std::unique_ptr<jtf::function::functionN1_t<double,int32_t> > lf(
            new jtf::function::least_square_Gauss_Newton(fptr));
      return lf.release();
    }

  }}
