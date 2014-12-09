#include "linear_solver.h"

#include <petscksp.h>
#include <petscversion.h>
#include <string>
#include <iostream>
#include <memory>
#include "./solver_pack/base/include/solver.h"
#include "dl_petsc_support.h"

using namespace hj::sparse;
//using namespace zjucad::matrix;
using namespace std;

class linear_solver_imp : public linear_solver
{
public:
	linear_solver_imp(hj::sparse::solver_base *slv)
		:slv_(slv) {
		assert(slv);
	}
	linear_solver_imp(PETsc_CG *p)
		:PETsc_CG_(p) {
	}

  virtual int solve(const double *b, double *x, size_t rhs, boost::property_tree::ptree &opts) {
		if(dynamic_cast<direct_solver_A *>(slv_.get()))
			slv_->solve(b, x, rhs);
		else if(PETsc_CG_.get()) {
			PETsc_CG_->solve(b, x, rhs, opts);
		}
		else {
		  opts.put("laspack/solver.desc","<CG,BiCGSTAB,...>");
		  opts.put("laspack/precond.desc","<Jacobi,SSOR,ILU>");
		  opts.put("laspack/iter.desc","number of iteration");
		  opts.put("laspack/precison.desc","precision");
		  opts.put("laspack/relax.desc","relax");
		  const string solve = opts.get("laspack/solver","CG"),
			precond = opts.get("laspack/precond", "SSOR");
			laspack_opts lopts(&solve[0], &precond[0]);
      lopts.iter_num_ = opts.get<size_t>("laspack/iter.value");
      lopts.precision_ = opts.get<double>("laspack/precison.value", 1e-5);
      lopts.relax_ = opts.get<double>("laspack/relax.value", 1.2);
			slv_->solve(b, x, rhs, &lopts);
		}
		return 0;
	}
protected:
	auto_ptr<hj::sparse::solver_base> slv_;
	auto_ptr<PETsc_CG> PETsc_CG_;
};

linear_solver *linear_solver::create(
					const double* val, const int32_t* idx,
					const int32_t * ptr, const size_t nnz,
					const size_t row, const size_t col,
				 	boost::property_tree::ptree & opts)
	{
		std::unique_ptr<hj::sparse::solver_base> slv;
		opts.put("linear_solver/type.desc","<PETsc, direct, iterative>");
		const string solver_type = opts.get<string>("linear_solver/type.value", "PETsc");
		opts.put("linear_solver/name.desc","<cholmod, umfpack>");

		if(solver_type == "direct")
			slv.reset(direct_solver_A::create(val, idx, ptr, nnz, row, col, 
																				opts.get<string>("linear_solver/name.value").c_str()));
		else if(solver_type == "iterative")
			slv.reset(iterative_solver_A::create(val, idx, ptr, nnz, row, col, "laspack"));
		else if(solver_type == "PETsc") {
			static auto_ptr<PETsc> petsc_init(ucreate_PETsc());
			return new linear_solver_imp(ucreate_PETsc_CG(val, idx, ptr, nnz, row, col, opts.get<string>("PETsc/pc.value", "sor").c_str()));
		}
		else
			cerr << "no such solver" << endl;
		if(slv.get()){
			return new linear_solver_imp(slv.release());
		}
		return 0;
	}
