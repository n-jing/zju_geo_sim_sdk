#include <string>
#include <iostream>

#include "func2opt.h"
#include "../include/optimizer.h"
#include "KKT.h"
#include "SQP.h"
#include "lbfgs.h"

using namespace std;

namespace jtf {

  int optimize(
      hj::math_func::math_func_t<double,int32_t> &f,
      zjucad::matrix::matrix<double> &x,
      boost::property_tree::ptree & pt,
      const std::vector<std::shared_ptr<hj::math_func::math_func_t<double,int32_t> > > *eqn_cons,
      const std::vector<std::shared_ptr<hj::math_func::math_func_t<double,int32_t> > > *ineqn_cons,
      jtf::opt::callbacks * cb){
  typedef shared_ptr<jtf::function::functionN1_t<double,int32_t> > jtf_func_ptr;
    func2opt<double,int32_t, RAW> fun(&f);
    unique_ptr<vector<jtf_func_ptr> > eqn_cons_jtf ;
    unique_ptr<vector<jtf_func_ptr> > ineqn_cons_jtf;
    if(eqn_cons){
        eqn_cons_jtf.reset(new vector<jtf_func_ptr>);
        for(size_t fi = 0; fi < eqn_cons->size(); ++fi){
            eqn_cons_jtf->push_back(jtf_func_ptr(new func2opt<double,int32_t,STD_SHARED>((*eqn_cons)[fi])));
          }
      }
    if(ineqn_cons){
        ineqn_cons_jtf.reset(new vector<jtf_func_ptr>);
        for(size_t fi = 0; fi < ineqn_cons->size(); ++fi){
            ineqn_cons_jtf->push_back(jtf_func_ptr(new func2opt<double,int32_t,STD_SHARED>((*ineqn_cons)[fi])));
          }
      }
    optimize(fun, x, pt, eqn_cons_jtf.get(), ineqn_cons_jtf.get(), cb);
    return 0;
  }

  int optimize(	jtf::function::functionN1_t<double,int32_t> &f,
                zjucad::matrix::matrix<double> &x,
                boost::property_tree::ptree &pt,
                const vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > * eqn_cons_,
                const vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > * ineqn_cons,
                jtf::opt::callbacks * cb)
  {
    using namespace std;
    pt.put("package.desc", "jtf");
    pt.put("alg.desc", "SQP/KKT");

    const string package = pt.get<string>("package.value");
    const string solver = pt.get<string>("alg.value");

    if(package != "jtf") {
      std::cerr << "# [error] package should be jtf." << std::endl;
      return __LINE__;
    }
    int rtn = 0;
    if(solver == "SQP"){
        if(eqn_cons_ !=  nullptr && eqn_cons_->size() > 1){
            cerr << "# [error] SQP do not support more than one equality constraints." << endl;
            return __LINE__;
          }
        SQP sqp;
        sqp.set_f(f);
        if(eqn_cons_){
            sqp.set_c(eqn_cons_->front().get());
          }
        if(cb)
          sqp.set_callbacks(cb);
        rtn = sqp.solve(&x[0], pt);
      }else if(solver == "SQPEX"){
        if(eqn_cons_ !=  nullptr && eqn_cons_->size() > 1){
            cerr << "# [error] SQPEX do not support more than one equality constraints." << endl;
            return __LINE__;
          }
        SQPEX sqpex;
        sqpex.set_f(f);
        if(eqn_cons_){
            sqpex.set_c(eqn_cons_->front().get());
          }
        rtn = sqpex.solve(&x[0], pt);
      }else if(solver == "KKT"){
        KKT k(f, eqn_cons_, ineqn_cons);
        rtn = k.solve(&x[0],pt);
      }else if(solver == "lbfgs"){
        cerr << "# [info] unconstraint optimization." << std::endl;
        lbfgs ls(f);
        rtn = ls.solve(&x[0], pt);
      } else {
        return __LINE__;
      }

    return rtn ;
  }
}
