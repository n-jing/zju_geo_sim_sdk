#include "./Gauss-Newton.h"

#include <time.h>
#include <string.h>

#include <iostream>
#include <fstream>

#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>

#include <hjlib/sparse/cached_AAT.h>
#include <zjucad/ptree/ptree.h>
//#include "../common/IO.h"
#include "fast_AAT.h"
#include "common.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;
using namespace hj::function;

#define DUMP_EQ 0

// double validate_jac(function &f, const double *x)
// {
// 	const itr_matrix<const double *> x0(f.dim_of_x(), 1, x);
// 	matrix<double> x1, f0(f.dim_of_f(), 1), f1 = f0;
// 	f.x_changed();
// 	f.val(&x0[0], &f0[0]);
	
// 	csc<double, int32_t> JT(f.dim_of_x(), f.dim_of_f(), f.jac_nnz()), J;
// 	f.jac(&x0[0], &JT.val()[0], &JT.ptr()[0], &JT.idx()[0]);
// 	trans(JT, J);

// 	const double eps = 1e-10;
// 	matrix<double> Ji;
// 	double rtn = 0;
// 	for(size_t xi = 0; xi < f.dim_of_x(); ++xi) {
// 		x1 = x0;
// 		x1[xi] += eps;
// 		f.x_changed();
// 		f.val(&x1[0], &f1[0]);
// 		Ji = (f1-f0)/eps;
// 		for(ptrdiff_t nzi = J.ptr()[xi]; nzi < J.ptr()[xi+1]; ++nzi)
// 			Ji[J.idx()[nzi]] -= J.val()[nzi];
// 		rtn += dot(Ji, Ji);
// 	}
// 	return rtn;
// }

Gauss_Newton::Gauss_Newton(const hj::function::function &f)
	:is_JT_sorted_(false)
{
	const ptrdiff_t nx = f.dim_of_x(), nf = f.dim_of_f();
	JT_.resize(f.dim_of_x(), f.dim_of_f(), f.jac_nnz());
	JT_.ptr()[0] = 0;
	g_.resize(nx);
	s_.resize(nx);
}

int Gauss_Newton::compute_g(const double *x, const hj::function::function &f, double *residual, zjucad::matrix::matrix<double> &g, double *norm2)
{
	itr_matrix<double *> f0(JT_.size(2), 1, residual);
	f.val(x, &f0[0], ctx_.get());
//	cerr << "compute_g eval jac from: " << x[0] << endl;
	f.jac(x, &JT_.val()[0], &JT_.ptr()[0], &JT_.idx()[0], ctx_.get());
#if DUMP_EQ
	{
		cerr << "dump JT, f beg" << endl;
		ofstream ofsJT("JT.csc", ofstream::binary);
		write_csc(ofsJT, JT_);
		ofstream ofsf("f.mat", ofstream::binary);
		matrix<double> f1 = f0;
		write_matrix(ofsf, f1);
		cerr << "dump JT, f end" << endl;
	}
#endif
	g(colon()) = 0;
	mv(false, JT_, f0, g);
	if(norm2) {
		norm2[0] = dot(f0, f0);
		norm2[1] = dot(g, g);
		cout << "residual: " << norm2[0]
			 << "\ngrad: " << norm2[1] << endl;
	}
	return 0;
}

int Gauss_Newton::solve_step(const double *x, double *s, boost::property_tree::ptree &opts)
{
	if(nnz(JTJ_))
		fast_AAT(JT_, JTJ_, is_JT_sorted_);
	else {
		is_JT_sorted_ = is_sorted_csc(JT_);
		cerr << "is_JT_sorted_: " << is_JT_sorted_ << endl;
		AAT<map_by_sorted_vector>(JT_, JTJ_);
	}

	slv_.reset(linear_solver::create(JTJ_, opts));
	if(!slv_.get()) return __LINE__;
	slv_->solve(&g_[0], s, 1, opts);
	itr_matrix<double *> s0(g_.size(), 1, s);
	s0 = -s0;
	return 0;
}

int Gauss_Newton::iterate(double *x, const hj::function::function &f, double *residual,boost::property_tree::ptree &opts)
{
	itr_matrix<double *> x0(s_.size(), 1, x);
	int iter_num = opts.get<int>("iter.value");
	for(int i = 0; i < iter_num; ++i) {
		compute_g(x, f, residual, g_);

		if(solve_step(x, &s_[0], opts))
			return __LINE__;

		x0 += s_;
		ctx_.reset(f.new_ctx(&x0[0]));
	}
	return 0;
}


damped_Gauss_Newton::damped_Gauss_Newton(const hj::function::function &f)
	:Gauss_Newton(f)
{
}

double simple_line_search(
	const double *x, const hj::function::function &f, const zjucad::matrix::matrix<double> &s, double prev_ratio)
{
	const itr_matrix<const double *> x0(s.size(), 1, x);
	matrix<double> x1(s.size()), f0(f.dim_of_f());
	auto_ptr<hj::function::func_ctx> ctx(f.new_ctx(&x0[0]));
	f.val(&x0[0], &f0[0], ctx.get());
	double residual[2] = {dot(f0, f0), 0};
	if(prev_ratio < 0.5)
		prev_ratio *= 2;
	for(;1;) {
		x1 = x0+prev_ratio*s;
		auto_ptr<hj::function::func_ctx> ctx1(f.new_ctx(&x1[0]));
		f.val(&x1[0], &f0[0], ctx1.get());
		residual[1] = dot(f0, f0);
		if(residual[1] > residual[0]) {
			if(prev_ratio < 1e-15) {
				cerr << "simple_line_search fail." << endl;
				return -1; // search fail
			}
			prev_ratio *= 0.5;
		}
		else
			break;
	}
	cout << "ratio: " << prev_ratio << endl;
	return prev_ratio;
}

int damped_Gauss_Newton::iterate(double *x, const hj::function::function &f, double *residual, boost::property_tree::ptree &opts)
{
	itr_matrix<double *> x0(s_.size(), 1, x);
	double ratio = 1.0;
	int iter_num = opts.get<int>("iter.value");
	for(int i = 0; i < iter_num; ++i) {
		compute_g(x, f, residual, g_);

		if(Gauss_Newton::solve_step(x, &s_[0], opts))
			return __LINE__;
		ratio = simple_line_search(x, f, s_, ratio);
		if(ratio < 0)
			return 1;
		x0 += ratio*s_;
	}
	return 0;
}

LM_Gauss_Newton::LM_Gauss_Newton(const hj::function::function &f)
	:Gauss_Newton(f)
{
}

template <typename VAL_TYPE, typename INT_TYPE>
inline void add_to_diag(csc<VAL_TYPE, INT_TYPE> &A, VAL_TYPE *diag)
{
	for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
		INT_TYPE nzi = A.ptr()[ci];
		for(; nzi < A.ptr()[ci+1]; ++nzi) {
			if(A.idx()[nzi] == ci) {
				A.val()[nzi] += diag[ci];
				break;
			}
		}
		if(nzi == A.ptr()[ci+1]) {
			cout << "bad add to diag." << endl;
		}
	}
}

template <typename VAL_TYPE, typename INT_TYPE>
inline void diag(const csc<VAL_TYPE, INT_TYPE> &A, matrix<VAL_TYPE> &d)
{
	d = zeros<VAL_TYPE>(A.size(1), 1);
	for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
		for(INT_TYPE nzi = A.ptr()[ci]; nzi < A.ptr()[ci+1]; ++nzi) {
			if(A.idx()[nzi] == ci)
				d[ci] = A.val()[nzi];
		}
	}
}

int LM_Gauss_Newton::solve_step(const double *x, double mu, double *s, boost::property_tree::ptree &opts)
{
	clock_t start;
	start = clock();
//	AAT<map_by_sorted_vector>(JT_, JTJ_);
//	static cached_AAT calc_AAT_;
//	calc_AAT_(JT_, JTJ_);
	if(nnz(JTJ_))
		fast_AAT(JT_, JTJ_, is_JT_sorted_);
	else {
		is_JT_sorted_ = is_sorted_csc(JT_);
		cerr << "is_JT_sorted_: " << is_JT_sorted_ << endl;
		AAT<map_by_sorted_vector>(JT_, JTJ_);
	}
	cerr << "AAT: " << clock()-start << endl;

	matrix<double> D(JTJ_.size(1));
	compute_D(&D[0]);
	D *= mu;
	add_to_diag(JTJ_, &D[0]);
	start = clock();
	slv_.reset(linear_solver::create(JTJ_, opts));
	cerr << "factor: " << clock()-start << endl;

	if(!slv_.get()) return __LINE__;
	start = clock();
	slv_->solve(&g_[0], s, 1, opts);
	cerr << "solve: " << clock()-start << endl;
	itr_matrix<double *> s0(g_.size(), 1, s);
	s0 = -s0;
	return 0;
}

int LM_Gauss_Newton::iterate(double *x, const hj::function::function &f, double *residual, boost::property_tree::ptree & opts)
{
	const double eps[3] = {
		opts.get<double>("epsg.value", 1e-8),
		opts.get<double>("epsf.value", 1e-8),
		opts.get<double>("epsx.value", 1e-8)
	};

	itr_matrix<double *> x0(s_.size(), 1, x), f0(JT_.size(2), 1, residual);
	double mu = 1e-3, alpha = 1;
	opts.put("GS_LM_mu.desc","mu for trust region.");
	if(zjucad::has("GS_LM_mu.value", opts))
		mu = opts.get<double>("GS_LM_mu.value");

	double norm2_of_f_g[2];

	compute_g(x, f, residual, g_, norm2_of_f_g);

	int iter_num = opts.get<int>("iter.value");
	opts.put("dump_step.desc","dump every n steps");
	size_t dump_step = opts.get<size_t>("dump_step.value", size_t(0));
	string dump_pref;
	opts.put("dump_pref.desc","prefix of the files for dumping x");
	if(dump_step > 0)
		dump_pref = opts.get<string>("dump_pref.value");
	for(int i = 0; i < iter_num; ++i) {
		if(norm2_of_f_g[1] < eps[0])
			break;
		const double current_residual = norm2_of_f_g[0];
		cout << "mu: " << mu << endl;
		if(solve_step(x, mu, &s_[0], opts)) {
      mu *= 4;
      continue;
    }

		if(1 || nnz(JT_) != f.jac_nnz()) { // in case of reduce memory, seems very fast
//			cerr << "beg recompute jac for reducing memory." << endl;
			JT_.resize(f.dim_of_x(), f.dim_of_f(), f.jac_nnz());
			JT_.ptr()[0] = 0;
//			cerr << "iterate eval jac from: " << x[0] << endl;
			f.jac(x, &JT_.val()[0], &JT_.ptr()[0], &JT_.idx()[0], ctx_.get());
//			cerr << trans(JT_.val()(colon(0, 8))) << endl;
//			cerr << "end recompute jac for reducing memory." << endl;
		}

		matrix<double> f1 = f0;
		mv(true, JT_, s_*alpha, f1);
		const double estimate_residual = dot(f1, f1);
		if(estimate_residual > current_residual) {
			cout << "ESTIMATE_RESIDUAL ERROR." << endl;
		}
		if(dot(s_, s_) < eps[2])
			break;
		x0 += s_;

		dump(i, dump_step, dump_pref.c_str(), x, x0.size());

		f.val(x, &f0[0], ctx_.get());
		const double real_residual = dot(f0, f0);
		double ratio = -1;
		if(real_residual < current_residual) {
			const double delta_f = real_residual-current_residual;
			if(delta_f < 0 && -delta_f < eps[1])
				break;
			ratio = delta_f/(estimate_residual-current_residual);
		}
		mu = adjust_mu(mu, ratio);
		compute_g(x, f, residual, g_, norm2_of_f_g);
	}
	return 0;
}

double LM_Gauss_Newton::adjust_mu(double mu, double ratio)
{
	const double min_max_mu[] = {1e-5, 1e6};
	if(ratio < 0) {
		mu *= 4;
		if(mu < min_max_mu[0])
			mu = min_max_mu[0];
		cout << "bad step." << endl;
	}
	else if(ratio < 0.25) {
		mu *= 4;
		if(mu < min_max_mu[0])
			mu = min_max_mu[0];
	}
	else
		mu /= sqrt(ratio*4);
	if(mu > min_max_mu[1])
		mu = min_max_mu[1];
	return mu;
}

void LM_Gauss_Newton::compute_D(double *D)
{
	matrix<double> diag_JTJ;
	diag(JTJ_, diag_JTJ);
	const double avg_diag = sum<double>(diag_JTJ)/diag_JTJ.size();
	for(ptrdiff_t i = 0; i < diag_JTJ.size(); ++i)
		D[i] = avg_diag;
}


More_Gauss_Newton::More_Gauss_Newton(const hj::function::function &f)
	:LM_Gauss_Newton(f)
{
	D_ = zeros<double>(f.dim_of_x(), 1);
}

int More_Gauss_Newton::solve_step(const double *x, double mu, double *s, boost::property_tree::ptree &opts)
{
	clock_t beg;
	beg = clock();
#if 0
//	AAT<map_by_sorted_vector>(JT_, JTJ_);
	static cached_AAT calc_AAT_; // cost too much memory with level 2 cache
	calc_AAT_(JT_, JTJ_, 1);
#else // fast_AAT is faster
	if(nnz(JTJ_))
		fast_AAT(JT_, JTJ_, is_JT_sorted_);
	else {
		is_JT_sorted_ = is_sorted_csc(JT_);
		cerr << "is sorted: " << is_JT_sorted_ << endl;
		AAT<map_by_sorted_vector>(JT_, JTJ_);
		cerr << JTJ_.size(1) << " " << JTJ_.size(2) << " " << nnz(JTJ_) << " "
			 << nnz(JTJ_)/double(JTJ_.size(1)) << " "
			 << nnz(JTJ_)/double(JTJ_.size(1)*JTJ_.size(2))*100 << "%" << endl;
	}
#endif
	cerr << "# JTJ time: " << double(clock()-beg)/CLOCKS_PER_SEC << endl;
#if DUMP_EQ
	{
		cerr << "begin dumping JTJ, JTf" << endl;
		ofstream ofsJTJ("JTJ.csc", ofstream::binary);
		write_csc(ofsJTJ, JTJ_);
		ofstream ofsJTf("JTf.mat", ofstream::binary);
		write_matrix(ofsJTf, g_);
		cerr << "end dumping JTJ, JTf" << endl;
	}
#endif

	compute_D(&D_[0]);
	{
		matrix<double> D = D_*mu;
		add_to_diag(JTJ_, &D[0]);
	}
//	cerr << trans(JT_.val()(colon(0, 8))) << endl;
	JT_.resize(0, 0, 0);

	itr_matrix<double *> s0(g_.size(), 1, s);
	beg = clock();
	slv_.reset(linear_solver::create(JTJ_, opts));
	cerr << "create time: " << double(clock()-beg)/CLOCKS_PER_SEC << endl;
	if(!slv_.get()) {
		cerr << "create linear solver fail." << endl;
		return __LINE__;
	}
	beg = clock();
	slv_->solve(&g_[0], s, 1, opts);
	cerr << "solve time: " << double(clock()-beg)/CLOCKS_PER_SEC << endl;
	slv_.reset(0);
#if 1
	matrix<double> Ax = zeros<double>(g_.size(), 1);
	mv(false, JTJ_, s, Ax);
	double r = norm(Ax-g_);
	cerr << "# |Ax-b|: " << r << " ||r||/||b||=" << r/norm(g_) << endl;
#endif
	s0 = -s0;
	return 0;
}

void More_Gauss_Newton::compute_D(double *D)
{
	matrix<double> D1 = zeros<double>(D_.size(), 1);
	for(size_t ci = 0; ci < JT_.size(2); ++ci) {
		for(size_t nzi = JT_.ptr()[ci]; nzi < JT_.ptr()[ci+1]; ++nzi) {
			D1[JT_.idx()[nzi]] += JT_.val()[nzi]*JT_.val()[nzi];
		}
	}
	D1 = sqrt(D1); // STRANGE: should be squared, DTD, according to
				   // the book, but removing the sqrt leads to slow
				   // convergence.
	for(ptrdiff_t i = 0; i < JT_.size(1); ++i) {
		if(D[i] < D1[i])
			D[i] = D1[i];
	}
}
