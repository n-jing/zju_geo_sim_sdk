#include "../src/optimization.h"

void function1_grad(const alglib::real_1d_array &x,

                    double &func,

                    alglib::real_1d_array &grad, void *ptr)

{

  //

  // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4

  // and its derivatives df/d0 and df/dx1

  //

  func = 100*pow(x[0]+3,4) + pow(x[1]-3,4);

  grad[0] = 400*pow(x[0]+3,3);

  grad[1] = 4*pow(x[1]-3,3);

}




int main()
{

  alglib::real_1d_array x = "[0,0]";

  double epsg = 0.0000000001;

  double epsf = 0;

  double epsx = 0;

  double stpmax = 0.1;

  alglib::ae_int_t maxits = 0;

  alglib::minlbfgsstate state;

  alglib::minlbfgsreport rep;



  // first run

  alglib::minlbfgscreate(1, x, state);

  alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);

  alglib::minlbfgssetstpmax(state, stpmax);

  alglib::minlbfgsoptimize(state, function1_grad);

  alglib::minlbfgsresults(state, x, rep);



  printf("%s\n", x.tostring(2).c_str()); // EXPECTED: [-3,3]



  // second run - algorithm is restarted

  x = "[10,10]";

  alglib::minlbfgsrestartfrom(state, x);

  alglib::minlbfgsoptimize(state, function1_grad);

  alglib::minlbfgsresults(state, x, rep);



  printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4

  printf("%s\n", x.tostring(2).c_str()); // EXPECTED: [-3,3]

  return 0;

}
