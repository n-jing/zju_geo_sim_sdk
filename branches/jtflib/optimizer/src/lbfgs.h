#ifndef JTF_LBFGS_H
#define JTF_LBFGS_H

#include "../include/opt.h"
#include "../../function/include/function.h"
#include <optimization.h>
#include <zjucad/ptree/ptree.h>

namespace jtf {
  class lbfgs
  {
  public:
    lbfgs(jtf::function::functionN1_t<double,int32_t> &func) : func_(func){}

    void eval_fg(const alglib::real_1d_array &x, double &func,
                 alglib::real_1d_array &grad);

    int solve(double* x, boost::property_tree::ptree &pt);

  private:
    jtf::function::functionN1_t<double,int32_t> &func_;
  };
}
#endif // LBFGS_H
