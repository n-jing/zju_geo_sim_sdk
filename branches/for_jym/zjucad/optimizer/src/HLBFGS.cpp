#include "HLBFGS.h"

#include "common.h"

#include <string>
#include <iostream>
#include <HLBFGS/HLBFGS.h>
#include <HLBFGS/Lite_Sparse_Matrix.h>
#include <zjucad/ptree/ptree.h>
#include <hjlib/sparse/sparse.h>

using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;
using namespace hj::sparse;

static const hj::function::function *g_fun = 0;
auto_ptr<matrix<double> > g_f;
auto_ptr<hj::sparse::csc<double, int32_t> > g_JT;
auto_ptr<Lite_Sparse_Matrix> g_H;

void HLBFGS_fg(int N, double* x, double *prev_x, double* f, double* g)
{
	auto_ptr<func_ctx> fc(g_fun->new_ctx(x));
	g_fun->val(x, &(*g_f)[0], fc.get());
	*f = dot(*g_f, *g_f);

	hj::sparse::csc<double, int32_t> &JT_ = *g_JT;

	g_fun->jac(x, &JT_.val()[0], &JT_.ptr()[0], &JT_.idx()[0], fc.get());
	double *pgrad = g;
	fill(pgrad, pgrad+g_fun->dim_of_x(), 0);
	mv(false, JT_, *g_f, pgrad);
}

void HLBFGS_h(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)
{
	static bool first = true;

	auto_ptr<func_ctx> fc(g_fun->new_ctx(x));

	if(first) {
		HLBFGS_fg(N, x, prev_x, f, g);
	}

	hj::sparse::csc<double, int32_t> &JT_ = *g_JT;
	g_fun->jac(x, &JT_.val()[0], &JT_.ptr()[0], &JT_.idx()[0], fc.get());

	hj::sparse::csc<double, int32_t> JTJ;
//	hj::sparse::AAT<map_by_sorted_vector>(JT_, JTJ);

	matrix<double> diag = zeros<double>(N, 1);
	for(size_t c = 0; c < JT_.size(2); ++c) {
		for(size_t nzi = JT_.ptr()[c]; nzi < JT_.ptr()[c+1]; ++nzi) {
			diag[JT_.idx()[nzi]] += JT_.val()[nzi]*JT_.val()[nzi];
		}
	}

	g_H.reset(new Lite_Sparse_Matrix(N, N,
									 SYM_LOWER, CCS, FORTRAN_TYPE, true));
	g_H->begin_fill_entry();

	// for(size_t ci = 0; ci < JTJ.size(2); ++ci) {
	// 	for(size_t nzi = JTJptr()[ci]; nzi < JTJptr()[ci+1]; ++nzi) {
	// 		if(JTJ.idx()[nzi] <= ci) continue;
	// 		g_H->fill_entry(JTJ.idx()[nzi], ci, JTJval()[nzi]);
	// 	}
	// }
	for(int i = 0; i < N; ++i)
	 	g_H->fill_diag(i, 1.0/sqrt(diag[i]));

	g_H->end_fill_entry();

	hessian.set_diag(g_H->get_diag());
	hessian.set_values(g_H->get_values());
	hessian.set_rowind(g_H->get_rowind());
	hessian.set_colptr(g_H->get_colptr());
	hessian.set_nonzeros(g_H->get_nonzero());

	first = false;

	cerr << "# set diag." << endl;
}

void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
	std::cout << iter <<":\t" << call_iter << "\t" << *f << "\t" << *gnorm  << std::endl;
}

int HLBFGS(const hj::function::function &fun, zjucad::matrix::matrix<double> &x,
		   zjucad::matrix::matrix<double> &residual,
		   boost::property_tree::ptree &opts)
{
	g_fun = &fun;
	g_f.reset(new matrix<double>(fun.dim_of_f(), 1));
	g_JT.reset(new hj::sparse::csc<double, int32_t>());
	g_JT->resize(fun.dim_of_x(), fun.dim_of_f(), fun.jac_nnz());

	int iter_num = opts.get<int>("iter.value");

	const int len = opts.get<int>("lbfgs-len.value");

	bool with_hessian = zjucad::has("hessian.value",opts);

	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	info[3] = 1;
	info[4] = iter_num;
	info[6] = with_hessian?10:0;
	info[7] = with_hessian?1:0;
	info[10] = 0;
	info[11] = 1;
 
	if(with_hessian) {
		g_H.reset(new Lite_Sparse_Matrix(fun.dim_of_x(), fun.dim_of_x(),
										 SYM_LOWER, CCS, C_TYPE, true));
		HLBFGS(fun.dim_of_x(), len,
			   &x[0], HLBFGS_fg, HLBFGS_h, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
		g_H.reset(0);
	}
	else {
		HLBFGS(fun.dim_of_x(), len,
			   &x[0], HLBFGS_fg, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
	}

	auto_ptr<func_ctx> fc(fun.new_ctx(&x[0]));
	fun.val(&x[0], &residual[0], fc.get());

	g_f.reset(0);
	g_JT.reset(0);

	return 0;
}
