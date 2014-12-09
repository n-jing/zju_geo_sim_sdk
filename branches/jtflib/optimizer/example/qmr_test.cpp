#include "../src/qmr.h"
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <iostream>
#include <functional>
#include "../../util/include/util.h"

using namespace std;
using namespace zjucad::matrix;

int test_qmr()
{
  double A[] = {1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
                0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1,
                0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, -1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
                0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0,
                0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, -1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
                1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, -1, -1, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0};
  double b[] = {-0, -0,  -0, -0, 2, -0, -0, -0, -0, 1,-1, 1, -1, -1, -1, -0,-0, -0, 0, 0,-2};

  matrix<double> A_m(21,21);
  matrix<double> b_m(21,1);

  std::copy(A,A+21*21, A_m.begin());
  std::copy(b,b+21, b_m.begin());

  matrix<double> x = zeros<double>(21,1);
  jtf::qmr_solver qs(A_m, b_m);
  qs.solve(&x[0], 100);
  cerr << norm(A_m * x - b_m) << endl;

  return 0;
}

int test_qmr2()
{
  matrix<double> A_m = eye<double>(10);
  A_m(9,9) = 0;
  matrix<double> b_m = ones<double>(10,1);

  matrix<double> x = zeros<double>(10,1);
  jtf::qmr_solver qs(A_m, b_m);
  qs.solve(&x[0], 100);

  cerr << norm(A_m * x - b_m) << endl;

  return 0;
}

class test {
public:
  test(){ p = new int[10];}
  test(const test & t){
    p = new int[10]; }
  test & operator = (const test &t){p = new int [10];}
 void print_info(double a)
  {
    cerr << a << endl;
    for(size_t i = 0; i < 1000000; ++i)
      a = cos(cos(a)/sin(a));
  }
  ~test(){
    cerr << __LINE__ << endl;
    delete []p;
  }
private:
  int *p;

};

// int test_qmr3()
// {
//   matrix<double> A_m = eye<double>(10);
//   matrix<double> b_m = rand<double>(10,1);

//   matrix<double> x = rand<double>(10,1);
//   hj::sparse::csc<double,int32_t> A_csc;

//   hj::sparse::convert(A_m, A_csc);

//   {
//     test t;
//     const std::function<void(void)> &f = std::bind(&test::print_info, std::cref(t), 2.0);
//     cout << jtf::util::func_time(f) << endl;
//   }

// //    jtf::PETsc_init pi;
// //    for(size_t i = 0; i < 2; ++i){
// //        jtf::PETsc_TFQMR qs(false, A_csc, b_m);
// //        std::function<int(void)> func = std::bind(&jtf::PETsc_TFQMR::solve, &qs,&x[0],50);
// //        cerr << "time " << jtf::util::func_time(func) << endl;
// //        //qs.solve(&x[0],50);
// //      }
// //    cerr << norm(A_m * x - b_m) << endl;

//   return 0;
// }
