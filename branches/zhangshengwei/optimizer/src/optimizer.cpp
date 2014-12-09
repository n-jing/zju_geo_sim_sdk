#include <string>
#include <iostream>

#include "../include/optimizer.h"
#include "KKT.h"
#include "SQP.h"
#include "lbfgs.h"

using namespace std;

namespace jtf {

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

    if(package != "jtf") return __LINE__;
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
