#include <jtflib/function/function.h>
#include <hjlib/sparse/sparse.h>
#include "opt.h"

namespace jtf{
  /// it solves min f(x)
  ///           Ax=b
  class KKT : public opt
  {
  public:
    KKT():f_(nullptr),eqn_cons_(nullptr), ineqn_cons_(nullptr){}
    KKT(jtf::function::functionN1_t<double,int32_t>  &f,
        const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > * eqn_cons,
        const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > * ineqn_cons)
      : f_(&f), eqn_cons_(eqn_cons), ineqn_cons_(ineqn_cons){}

    virtual ~KKT(){}

    virtual int set_f(jtf::function::functionN1_t<double,int32_t> &f){
      f_ = &f;
    }

    virtual int set_c(const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > * eqn_cons,
                      const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > * ineqn_cons)
    {
      eqn_cons_ = eqn_cons;
      ineqn_cons_ = ineqn_cons;
      return 0;
    }
    virtual int solve(double *x, boost::property_tree::ptree &pt);

    int get_A_b(hj::sparse::csc<double,int32_t> &A,
                zjucad::matrix::matrix<double> &b) const;
  protected:
//    int preprocess_init_data(double *x,boost::property_tree::ptree &pt);

    int solve_sparse(double *x, boost::property_tree::ptree &pt);

    int solve_alglib(double *x, boost::property_tree::ptree &pt);

  private:
    jtf::function::functionN1_t<double,int32_t> * f_;
    std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > combined_obj_;
    const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > *eqn_cons_;
    const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > *ineqn_cons_;
  };
}
