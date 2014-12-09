#include <iostream>
#include <linear_solver.h>

using namespace std;

int main(int argc, char *argv[])
{
  /**
   * a=
   1   2   0   0
   0   2   0   0
   0   0   2   0
   0   0   0   1
   x =
   b =
   4
   4
   3
   1
   * 
   */
   cout << "into main" << endl;
  zjucad::matrix::matrix<double> mat = zjucad::matrix::rand<double>(4,4);
  cout << "init finish" << endl;
  // // mat(0,0) = 1; mat(0,1) = 2; mat(0,2) = 0; mat(0,3) = 0;
  // // mat(1,0) = 0; mat(1,1) = 2; mat(1,2) = 0; mat(1,3) = 0;
  // // mat(2,0) = 1; mat(2,1) = 0; mat(2,2) = 2; mat(2,3) = 0;
  // // mat(3,0) = 0; mat(3,1) = 0; mat(3,2) = 0; mat(3,3) = 1;
  mat = temp(trans(temp(mat)) * mat);

  for (int i=0; i<4; ++i)
  {
    for (int j=0; j<4; ++j)
    {
      cout << mat(i,j) << " ";
    }
    cout << endl;
  }
  hj::sparse::csc<double, int32_t> csc_mat;
  boost::property_tree::ptree opts;
  // opts.put("linear_solver/type.value", "PETsc");
  opts.put("linear_solver/type.value", "direct");
  // opts.put("linear_solver/name.value", "cholmod");
  opts.put("linear_solver/name.value", "umfpack");
  convert(mat, csc_mat, 1e-4);
  linear_solver* solver = linear_solver::create(csc_mat, opts);
  zjucad::matrix::matrix<double> b(4,1);
  b[0] = 4;
  b[1] = 4;
  b[2] = 3;
  b[1] = 1;
  zjucad::matrix::matrix<double> x(4,1);

  boost::property_tree::ptree opts_solve;
  solver->solve(&b[0],&x[0],1,opts_solve);

  zjucad::matrix::matrix<double> diff = -1*b;
  hj::sparse::mv(false, csc_mat, x, diff);
  cerr << norm(diff) << endl;
  
  delete solver;
  return 0;
}

