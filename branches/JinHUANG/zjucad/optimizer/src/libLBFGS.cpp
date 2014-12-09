#include "libLBFGS.h"

#include "common.h"

#include <limits>

extern "C" {
#include <libLBFGS/lbfgs.h>
}

#include <iostream>
#include <hjlib/sparse/sparse.h>

//using namespace hj::function;
using namespace hj::sparse;
using namespace std;

class libLBFGS_ctx
{
public:
	libLBFGS_ctx(const hj::function::function &f, size_t nrow, size_t ncol = 1)
		:fun_(f), current_iter_(0), dump_step_(0), nrow_(nrow), ncol_(ncol) {
		f_.resize(fun_.dim_of_f());
		JT_.resize(fun_.dim_of_x(), fun_.dim_of_f(), fun_.jac_nnz());
	}
	const hj::function::function &fun_;
	zjucad::matrix::matrix<double> f_;
	hj::sparse::csc<double, int32_t> JT_;
	size_t current_iter_;
	const size_t nrow_, ncol_;
	size_t dump_step_;
	string dump_pref_;
};

static lbfgsfloatval_t libLBFGS_evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
	libLBFGS_ctx *ctx = reinterpret_cast<libLBFGS_ctx *>(instance);
	auto_ptr<hj::function::func_ctx> fc(ctx->fun_.new_ctx(x));
	ctx->fun_.val(x, &ctx->f_[0], fc.get());

	ctx->fun_.jac(x, &ctx->JT_.val()[0], &ctx->JT_.ptr()[0], &ctx->JT_.idx()[0], fc.get());

	fill(g, g+ctx->fun_.dim_of_x(), 0);
	mv(false, ctx->JT_, ctx->f_, g);

	return dot(ctx->f_, ctx->f_);	
}

static int libLBFGS_progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
	libLBFGS_ctx *ctx = reinterpret_cast<libLBFGS_ctx *>(instance);
	dump(ctx->current_iter_, ctx->dump_step_, ctx->dump_pref_.c_str(),
		 x,
		 ctx->fun_.dim_of_x());
	++ctx->current_iter_;

    cerr << "# Iteration: " << k
		 << ", fx = " << fx
		 << ", xnorm = " << xnorm
		 << ", gnorm = " << gnorm
		 << ", step = " << step << "\n";

    return 0;
}

int libLBFGS(const hj::function::function &fun, zjucad::matrix::matrix<double> &x,
			 zjucad::matrix::matrix<double> &residual,
			 boost::property_tree::ptree &opts)
{
	libLBFGS_ctx ctx(fun, x.size(1), x.size(2));
	opts.put("dump_step.desc","dump step");
	ctx.dump_step_ = opts.get<size_t>("dump_step", size_t(0));
	opts.put("dump_pred,desc","dump pref");
	if(ctx.dump_step_ > 0)
	  ctx.dump_pref_ = opts.get<string>("dump_pref");

	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);
//	param.xtol = numeric_limits<double>::epsilon();
	param.epsilon = opts.get<double>("libLBFGS/eps.value", 1e-8);
	param.max_iterations = opts.get<int>("iter.value");

	int ret = lbfgs(fun.dim_of_x(), &x[0], 0,
					libLBFGS_evaluate, libLBFGS_progress,
					&ctx, &param);
	auto_ptr<hj::function::func_ctx> fc(ctx.fun_.new_ctx(&x[0]));
	ctx.fun_.val(&x[0], &residual[0], fc.get());

	return ret;
}

