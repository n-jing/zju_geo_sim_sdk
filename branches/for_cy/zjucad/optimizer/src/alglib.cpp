#include "alglib.h"

#include "common.h"

#include <string>
#include <iostream>

#include <optimization.h>

#include <hjlib/sparse/sparse.h>
#include <zjucad/ptree/ptree.h>
using namespace alglib;
using namespace zjucad::matrix;
using namespace std;

class alglib_ctx
{
public:
	alglib_ctx(const hj::function::function &f, size_t nrow, size_t ncol = 1)
		:fun_(f), current_iter_(0), dump_step_(0), nrow_(nrow), ncol_(ncol) {
		f_.resize(fun_.dim_of_f());
		JT_.resize(fun_.dim_of_x(), fun_.dim_of_f(), fun_.jac_nnz());
	}

	virtual ~alglib_ctx(){}

	virtual void eval_fg(const real_1d_array &x, double &func, real_1d_array &grad) {
		auto_ptr<hj::function::func_ctx> fc(this->fun_.new_ctx(x.getcontent()));
		this->fun_.val(x.getcontent(), &this->f_[0], fc.get());
		func = dot(this->f_, this->f_);

		this->fun_.jac(x.getcontent(), &this->JT_.val()[0], &this->JT_.ptr()[0], &this->JT_.idx()[0], fc.get());

		double *pgrad = grad.getcontent();
		fill(pgrad, pgrad+this->fun_.dim_of_x(), 0);
		mv(false, this->JT_, this->f_, pgrad);
	}

	virtual void rep(const real_1d_array &x, double func) {
		dump(this->current_iter_, this->dump_step_, this->dump_pref_.c_str(),
			 x.getcontent(),
			 this->fun_.dim_of_x());
    cout << "residual: " << func << endl;
		++this->current_iter_;
	}

	const hj::function::function &fun_;
	zjucad::matrix::matrix<double> f_;
	hj::sparse::csc<double, int32_t> JT_;
	zjucad::matrix::matrix<double> diag_H_;
	size_t current_iter_;
	const size_t nrow_, ncol_;
	size_t dump_step_;
	string dump_pref_;
};

class alglib_lbfgs_ctx : public alglib_ctx
{
public:
	alglib_lbfgs_ctx(const hj::function::function &f, size_t nrow, size_t ncol = 1)
		:alglib_ctx(f, nrow, ncol), state_(0) {
	}
	virtual void eval_fg(const real_1d_array &x, double &func, real_1d_array &grad) {
		alglib_ctx::eval_fg(x, func, grad);
		if(r1a_diag_H_.get() && state_) {
			diag_H_(colon()) = 0;
			for(size_t c = 0; c < JT_.size(2); ++c) {
				for(size_t nzi = JT_.ptr()[c]; nzi < JT_.ptr()[c+1]; ++nzi) {
					diag_H_[JT_.idx()[nzi]] += JT_.val()[nzi]*JT_.val()[nzi];
				}
			}
			static const double eps = 1e-8;
			for(size_t i = 0; i < diag_H_.size(); ++i) {
				if(diag_H_[i] < eps) {
					cerr << "# 0 in diag_H_." << endl;
					diag_H_[i] = eps;
				}
			}
			r1a_diag_H_->setcontent(diag_H_.size(), &diag_H_[0]);
		}
	}
	virtual void rep(const real_1d_array &x, double func) {
		alglib_ctx::rep(x, func);
		if(state_)
			minlbfgssetprecdiag(*state_, *r1a_diag_H_);
	}
	minlbfgsstate *state_;
	auto_ptr<real_1d_array> r1a_diag_H_;
};

void function2alglib_fg(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
	alglib_ctx *ctx = reinterpret_cast<alglib_ctx *>(ptr);
	ctx->eval_fg(x, func, grad);
}

void  alglib_rep(const real_1d_array &x, double func, void *ptr)
{
	alglib_ctx *ctx = reinterpret_cast<alglib_ctx *>(ptr);
	ctx->rep(x, func);
}

int alglib_min(const hj::function::function &fun, zjucad::matrix::matrix<double> &x,
			   zjucad::matrix::matrix<double> &residual,
			   boost::property_tree::ptree &opts)
{
	int iter_num = opts.get<int>("iter.value");
	if(iter_num <= 0)
		return 0;

	real_1d_array alglib_x;
	alglib_x.setcontent(x.size(), &x[0]);

	int rtn;

	const double eps[3] = {
		opts.get<double>("epsg.value", 1e-8),
		opts.get<double>("epsf.value", 1e-8),
		opts.get<double>("epsx.value", 1e-8)
	};

	auto_ptr<alglib_ctx> ctx;//(new alglib_ctx(fun, x.size(1), x.size(2)));
	if(opts.get<string>("alg.value") == "non-linear-cg") {
		ctx.reset(new alglib_ctx(fun, x.size(1), x.size(2)));
	}
	else if(opts.get<string>("alg.value") == "lbfgs") {
		ctx.reset(new alglib_lbfgs_ctx(fun, x.size(1), x.size(2)));
	}
	opts.put("dump_step.desc","dump step");
	ctx->dump_step_ = opts.get<size_t>("dump_step.value", size_t(0));
	if(ctx->dump_step_ > 0)
	  {
		opts.put("dump_pref.desc","dump pref");
		ctx->dump_pref_ = opts.get<string>("dump_pref.value");
	  }

	if(opts.get<string>("alg.value") == "non-linear-cg") {
		mincgstate state;
		mincgreport rep;

		mincgcreate(alglib_x, state);
		mincgsetxrep(state, true);
		mincgsetcond(state, eps[2], eps[1], eps[0], iter_num);
		mincgsetcgtype(state, -1);
		mincgoptimize(state, function2alglib_fg, alglib_rep, ctx.get());
		mincgresults(state, alglib_x, rep);

		rtn = int(rep.terminationtype);
	}
	else if(opts.get<string>("alg.value") == "lbfgs") {
		minlbfgsstate state;
		minlbfgsreport rep;

		if(zjucad::has("lbfgs-precdiag.value",opts)) {
      cerr << "# use precdiag" << endl;
			ctx->diag_H_.resize(x.size(), 1);
			alglib_lbfgs_ctx *p = dynamic_cast<alglib_lbfgs_ctx *>(ctx.get());
			if(!p) {
				cerr << "# bad alglib_ctx." << endl;
				return __LINE__;
			}
			p->state_ = &state;
			p->r1a_diag_H_.reset(new real_1d_array);
		}

		const int len = opts.get<int>("lbfgs-len.value");
		minlbfgscreate(x.size(), len, alglib_x, state);
		minlbfgssetxrep(state, true);
		minlbfgssetcond(state, eps[2], eps[1], eps[0], iter_num);
//		minlbfgssetdefaultpreconditioner(state); // default has been set, so this function can be ignored
		try {
			minlbfgsoptimize(state, function2alglib_fg, alglib_rep, ctx.get());
		}
		catch (ap_error ae) {
      cerr << "# ap_error: " << ae.msg << endl;
			return __LINE__;
		}
		minlbfgsresults(state, alglib_x, rep);
		rtn = int(rep.terminationtype);
	}
	// else if(!strcmp(solver_name, "minlm")) { // seems it only handle dense H
	// 	return __LINE__;
	// }
	else {
    cerr << "no such solver in alglib." << endl;
		return __LINE__;
	}
	copy(alglib_x.getcontent(), alglib_x.getcontent()+x.size(), &x[0]);
	auto_ptr<hj::function::func_ctx> fc(ctx->fun_.new_ctx(&x[0]));
	ctx->fun_.val(&x[0], &residual[0], fc.get());

	return rtn;
}
