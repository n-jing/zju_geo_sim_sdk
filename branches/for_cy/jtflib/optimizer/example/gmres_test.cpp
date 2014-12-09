#include "../src/gmres.h"
#include "../src/precondition.h"
#include <iostream>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>

using namespace std;
using namespace zjucad::matrix;

int test_gmres1()
{
  matrix<double> A = rand<double>(10,10);
  matrix<double> b = rand<double>(10,1);
  matrix<double> x = rand(10,1);

//  cerr << "A = " << A << endl;
//  cerr << "b = " << b << endl;

  jtf::mesh::read_matrix("A.mat", A);
  jtf::mesh::read_matrix("x.mat", x);
  jtf::mesh::read_matrix("b.mat", b);

//    cerr << "x0 = " << x << endl;
//    x[0] =   0.581071;
//    x[1] =-0.0195524;
//    x[2] =-0.286983;
//    x[3] =-0.525714;
//    x[4] =-0.261275;
//    x[5] =-0.121119;
//    x[6] = 0.54472;
//    x[7] =0.465294;
//    x[8] =-0.912155;
//    x[9] =1.01387;

//  jtf::mesh::write_matrix("A.mat", A);
//  jtf::mesh::write_matrix("x.mat", x);
//  jtf::mesh::write_matrix("b.mat", b);

  jtf::gmres_solver gs(A,b);

  std::shared_ptr<jtf::jacobi_precondition> jp(new jtf::jacobi_precondition(A,b));
  for(size_t i = 0; i < 5; ++i)
  gs.solve2(&x[0],5,5, jp);

  //gs.solve(&x[0],10,10);

  //cerr << x << endl;
  cerr << norm(A * x - b) << endl;
  cerr << b.size() << endl;
  return 0;
}

int test_gmres2()
{
  size_t N = 10;
  const double eps = 0.05;
  matrix<double> A = eye<double>(N) * 2.1;
  for(size_t i = 0; i < N; ++i){
      if(i > 0) A(i-1,i) = -1 + eps;
      if(i+1<A.size(1)) A(i+1,i) = -1-eps;
    }

  matrix<double> exact = ones<double>(N,1);

  matrix<double> b = A * exact;
  matrix<double> x = b;

  cerr << "A = " << A << endl;
  cerr << "b = " << b << endl;
  cerr << "x = " << x << endl;

  jtf::gmres_solver gs(A,b);

  std::shared_ptr<jtf::jacobi_precondition> jp(new jtf::jacobi_precondition(A,b));

  gs.solve2(&x[0],5,35, jp);

  //for(size_t i = 0; i < 7; ++i)
  //gs.solve(&x[0],5,35);

  //cerr << x << endl;
  cerr << norm(A * x - b) << endl;
  cerr << b.size() << endl;
  return 0;
}


int test_gmres3()
{

	hj::sparse::csc<double,int32_t> A(3,3,4);
	
	A.ptr()[1]	 = A.ptr()[0] + 2;
	A.idx()[A.ptr()[0] + 0] = 0;
	A.val()[A.ptr()[0] + 0] = 1;

	A.idx()[A.ptr()[0] + 1] = 2;
	A.val()[A.ptr()[0] + 1] = 1;

	A.ptr()[2]	 = A.ptr()[1] + 1;
	A.idx()[A.ptr()[1] + 0] = 1;
	A.val()[A.ptr()[1] + 0] = 1;

	A.ptr()[3]	 = A.ptr()[2] + 1;
	A.idx()[A.ptr()[2] + 0] = 2;
	A.val()[A.ptr()[2] + 0] = 1;

	matrix<double> b = rand<double>(3,1);

	std::shared_ptr<jtf::jacobi_precondition_sparse> jp(new jtf::jacobi_precondition_sparse(A,b));
	jtf::gmres_solver_sparse gs(A,b);

	matrix<double> r = -1*b;
	matrix<double> x = b*2;
	//	gs.solve2(&x[0], 1,10, jp);
	gs.solve(&x[0], 1,10);
	
	hj::sparse::mv(false, A, x, r);
	cerr << x << endl;
	cerr << "# [info] Ax-b " << norm(r) << endl;
	return 0;
}


