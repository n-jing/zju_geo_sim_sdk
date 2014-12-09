#ifndef JTF_SQP_H_
#define JTF_SQP_H_

#include "../include/opt.h"
#include <zjucad/optimizer/linear_solver.h>

namespace jtf {

  class SQP : public opt
  {
  public:
    SQP();
    virtual int set_f(jtf::function::functionN1_t<double,int32_t> &f);
    virtual int set_c(jtf::function::functionN1_t<double,int32_t> *c);
    virtual int solve(double *x, boost::property_tree::ptree &pt);
  protected:
    int solve_by_xx_decomposition(double *x, boost::property_tree::ptree &pt);
    jtf::function::functionN1_t<double,int32_t> *f_, *c_;
    std::unique_ptr<linear_solver> slv_;
    hj::sparse::csc<double, int32_t> H_, corrected_H_;
    zjucad::matrix::matrix<double> g_, s_, D_, tmp_;
  };

  class SQPEX:public SQP
  {
  public:
      SQPEX();
      virtual int solve(double *x, boost::property_tree::ptree &pt);
  protected:
      int solve_by_xx_decomposition(double *x, boost::property_tree::ptree &pt);
  };

}

#endif
