#include "lbfgs.h"
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/matrix.h>

namespace jtf {

  void function2alglib_fg(const alglib::real_1d_array & x, double &func,
                          alglib::real_1d_array &grad, void *ptr)
  {
    lbfgs *ls = reinterpret_cast<lbfgs*>(ptr);
    ls->eval_fg(x,func,grad);
  }

  void lbfgs::eval_fg(const alglib::real_1d_array &x, double &func,
                      alglib::real_1d_array &grad)
  {
    func = 0;
    func_.val(&x[0], func);
    zjucad::matrix::itr_matrix<double*> gra(func_.dim(),1,grad.getcontent());
    gra *= 0;
    func_.gra(&x[0], &gra[0]);
    std::cerr << "resid: " << func << "\t gra_norm " << zjucad::matrix::norm(gra) << std::endl;
  }

  int lbfgs::solve(double *x, boost::property_tree::ptree &pt)
  {
    int iter_num = pt.get<int>("iter.value");
    if(iter_num <= 0) return 0;

    alglib::real_1d_array alglib_x;
    alglib_x.setcontent(func_.dim(),x);

    const double eps[3] = {
      pt.get<double>("epsg.value",1e-12),
      pt.get<double>("epsf.value",1e-12),
      pt.get<double>("epsx.value",1e-12)
    };

    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;

    const int len = pt.get<int>("lbfgs-len.value");
    alglib::minlbfgscreate(func_.dim(), len, alglib_x, state);
    alglib::minlbfgssetxrep(state, true);
    alglib::minlbfgssetcond(state, eps[2], eps[1], eps[0], iter_num);
    try{
      alglib::minlbfgsoptimize(state, function2alglib_fg, nullptr, this);
    }catch(alglib::ap_error &e){
      std::cerr << e.msg << std::endl;
    }

    alglib::minlbfgsresults(state, alglib_x, rep);
    int rtn = int(rep.terminationtype);
    std::copy(alglib_x.getcontent(), alglib_x.getcontent()+func_.dim(), x);

    return rtn;
  }
}
