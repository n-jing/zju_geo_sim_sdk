#include "../include/function.h"
#include "neg_log.h"
#include <hjlib/function/function.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>

namespace jtf{
  namespace function{

    jtf::function::functionN1_t<double,int32_t> *
    neg_log_warpper(jtf::function::functionN1_t<double,int32_t> &f)
    {
      std::unique_ptr<jtf::function::functionN1_t<double,int32_t> > p(
            new neg_log_func(f));
      return p.release();
    }

    jtf::function::functionN1_t<double,int32_t> *
    neg_log_warpper(std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > fptr)
    {
      assert(fptr);
      std::unique_ptr<jtf::function::functionN1_t<double,int32_t> > p(
            new neg_log_func(fptr));
      return p.release();
    }

    jtf::function::functionN1_t<double,int32_t> *
    neg_log_warpper(std::shared_ptr<hj::function::function_t<double,int32_t> > hptr,
                    const zjucad::matrix::matrix<double> & weight)
    {
      std::unique_ptr<jtf::function::functionN1_t<double,int32_t> > p(
            new neg_log_from_hj_func(hptr, weight));
      return p.release();
    }

  }
}
