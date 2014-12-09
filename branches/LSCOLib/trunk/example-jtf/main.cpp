//#include "kernel.h"
//#include "direct_solver_kernel.h"

#include "kernel_common_file_pzr.h"
#include "direct_solver_kernel_common_file_pzr.h"

#include "kernel_zjucad_matrix.h"
#include "direct_solver_kernel_zjucad_matrix.h"

//#include "kernel_debug.h"
//#include "quasi_newton_debug.h"
//#include "problem_largest_small_polygon.h"
//#include "problem_smooth.h"
//#include "problem_simple_function.h"
#include "problem_arap.h"
//#include "problem_frame_align.h"
using namespace zjucad::LSCO;

int main(int argc, char * argv[])
{
    //    //kernel_debug<double>();

    //    //typedef kernel_traits<plain_matrix<double> > kernel_type;
    //    //typedef kernel_traits<COMMON::FixedSparseMatrix<double,COMMON::Kernel<double> > > kernel_type;
    //    typedef kernel_traits<hj::sparse::csc<double> > kernel_type;
    //    //quasi_newton_debug<kernel_type>(true,true);
    //test_smoothing_hard_constraint(argv[1]);
    //hanging_chain(1.0,3.0,4,200);
    // tet_smoothing(argv[1]);
    test_arap(argv[1], argv[2], argv[3]);
    //test_smoothing_soft_constraint(argv[1]);
    //test_simple();
    //largest_small_polygon(4);
    //test_align_frame(argv[1], argv[2], argv[3]);
    //    //system("pause");
}

/*
#include "kernel_common_file_pzr.h"
#include "direct_solver_kernel_common_file_pzr.h"
#include "trust_region_solver_steihuag_cg.h"
#include "kernel_debug.h"
#include <fstream>

using namespace zjucad::LSCO;

int main()
{
	typedef kernel_traits<COMMON::FixedSparseMatrix<double,COMMON::Kernel<double> > > kernel_type;
	//typedef kernel_traits<plain_matrix<double> > kernel_type;
	typedef kernel_type::value_type value_type;
	typedef kernel_type::vector_type vector_type;
	typedef kernel_type::sparse_matrix_type sparse_matrix_type;

	boost::property_tree::ptree pt;
	pt.put<bool>("trust_region_solver_in.report",true);
	pt.put<bool>("trust_region_solver_in.nonsingular_preconditioner",true);
	trust_region_solver_steihaug_cg<kernel_type> tr(pt);

	boost::shared_ptr<sparse_matrix_type> A(new sparse_matrix_type);
	boost::shared_ptr<sparse_matrix_type> s(new sparse_matrix_type);
	boost::shared_ptr<sparse_matrix_type> M(new sparse_matrix_type);
	boost::shared_ptr<sparse_matrix_type> invM(new sparse_matrix_type);
	vector_type g;
	vector_type x;
	A->resize(1,20);
	s->resize(20,20);
	M->resize(20,20);
	invM->resize(20,20);
	x.resize(20);
	g.resize(20);
	for(int i=0;i<20;i++)
	{
		s->addToElement(i,i,rand()*2.0f/(value_type)RAND_MAX-0.2f);
		M->addToElement(i,i,sqrt((*s)(i,i)));
		invM->addToElement(i,i,1.0f/sqrt((*s)(i,i)));
		g[i]=rand()/(value_type)RAND_MAX;
	}
	A->addToElement(0,0,1.0f);
	A->addToElement(0,5,-5.0f);

	boost::shared_ptr<krylov_matrix<kernel_type> > sk(new default_krylov_matrix<kernel_type>(*s));
	tr.set_hessian(sk);

	tr.set_constraint(A);

	//tr.set_preconditioner(M);
	//tr.set_preconditioner_inv(invM);
	//boost::shared_ptr<krylov_matrix<kernel_type> > invMk(new default_krylov_matrix<kernel_type>(*invM));
	//tr.set_preconditioner_inv(invMk);

	tr.assemble();
	tr.solve(g,100.0f,x);
	std::cout << x[0]+x[5]*-5.0f;
	system("pause");
}

#include "kernel.h"
#include "direct_solver_kernel.h"
#include "objective_function.h"
#include <boost/shared_ptr.hpp>

using namespace zjucad::LSCO;

int main()
{
	//kernel_debug<double>();
	typedef kernel_traits<plain_matrix<double> > kernel_type;
	boost::shared_ptr<objective_function<kernel_type> > kernel;
	//kernel.reset(new objective_function_L_D_BFGS<kernel_type>());
	kernel.reset(new objective_function_L_SR1<kernel_type>());
	system("pause");
}
}
*/

//#include "LSCO_solver_interior_cg.h"

//using namespace zjucad::LSCO;

//int main()
//{
//    typedef kernel_traits<plain_matrix<double> > kernel_type;
//    boost::property_tree::ptree pt;
//    LSCO_solver_interior_cg<kernel_type> sol(pt);

//    objective_function<kernel_type> obj;
//    kernel_type::vector_type x;
//    sol.solve(obj,x);
//}
