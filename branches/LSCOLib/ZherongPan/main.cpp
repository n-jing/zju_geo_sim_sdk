#include "kernel_common_file_pzr.h"
#include "direct_solver_kernel_common_file_pzr.h"
#include "LSCO_solver_interior_cg.h"
#include "LSCOProblemQP.h"

using namespace zjucad::LSCO;

int main(int argc, char * argv[])
{
#define NR 50

	ProblemQP qp(NR,0,10,false);

	boost::property_tree::ptree pt;
	pt.put<bool>("LSCO_solver_in.report",true);
	LSCO_solver_interior_cg<PROBLEM_KERNEL> sol(pt);
	
	PROBLEM_KERNEL::vector_type x;
	x.resize(NR);
	x.setZero();
	//x=qp._x;
	sol.solve(qp,x);
	system("pause");
}

/*#include "kernel_common_file_pzr.h"
#include "direct_solver_kernel_common_file_pzr.h"
#include "trust_region_solver_debug.h"

using namespace COMMON;
using namespace zjucad::LSCO;

int main()
{
	typedef kernel_traits<FixedSparseMatrix<scalarD,Kernel<scalarD> > > kernel_type;
	trust_region_solver_debug<kernel_type>(true,true,false,true);
	system("pause");
}*/

/*#include <CommonFilePZR/solvers/LinearSolver.h>

using namespace COMMON;

int main()
{
	FixedSparseMatrix<scalarD,Kernel<scalarD> > lhs;
	Kernel<scalarD>::Vec rhs,x;

	lhs.read(ifstream("lhs.dat",ios::binary));
	readBinaryData(rhs,ifstream("rhs.dat",ios::binary));

	boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > > cb(new Callback<scalarD,Kernel<scalarD> >);
	PCGSolver<scalarD,Kernel<scalarD>,NoPreconSolver<scalarD> > sol;
	sol.setCallback(cb);
	sol.setMatrix(lhs,true);
	sol.setSolverParameters(1E-6f,10000);
	sol.solve(rhs,x);
	system("pause");
}*/