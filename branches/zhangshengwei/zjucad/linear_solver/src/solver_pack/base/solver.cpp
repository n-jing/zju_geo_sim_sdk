#include <string.h>
#include <inttypes.h>

#include <string>
#include <memory>
#include <iostream>

#include <hjlib/sparse/sparse.h>
#include "solver_ctx.h"
#include "include/solver.h"
#include "../laspack/laspack.h"

using namespace std;

void *hj_sparse_direct_solver_A_create(
	unsigned char sizeof_int, 	unsigned char sizeof_val, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	const void *ptr, const void *idx, const void *val,
	const char *name, void *opts)
{
	if(!name)
		name = "umfpack";
	if(!strcmp(name, "umfpack")) {
		if(sizeof_val != sizeof(double))// umfpack only support double precision
			return 0;
		auto_ptr<umfpack_ctx> uc(
			new umfpack_ctx(
				sizeof_int, real_or_complex,
				nrows, ncols,
				ptr, idx, reinterpret_cast<const double *>(val), opts)
			);
		if(uc.get() && !uc->init())
			return uc.release();
	}
            else if(!strcmp(name, "cholmod")) { //for compatibility consideration, we name the default solver 'ldlt' method as 'cholmod'
		if(sizeof_val != sizeof(double))// cholmod only support double precision
			return 0;
		auto_ptr<cholmod_ctx> c(
			new cholmod_ctx(
				sizeof_int, sizeof_val, real_or_complex,
				nrows, ncols,
				ptr, idx, reinterpret_cast<const double *>(val), opts)
            );
        if(c.get() && !c->init(CHOLMOD_SIMPLICIAL))
			return c.release();
    }else if(!strcmp(name, "cholmod_llt")) {
        if(sizeof_val != sizeof(double))// cholmod only support double precision
            return 0;
        auto_ptr<cholmod_ctx> c(
            new cholmod_ctx(
                sizeof_int, sizeof_val, real_or_complex,
                nrows, ncols,
                ptr, idx, reinterpret_cast<const double *>(val), opts)
            );
        if(c.get() && !c->init(CHOLMOD_SUPERNODAL)){
            return c.release();
        }
    }else if(!strcmp(name, "cholmod_auto")) {
        if(sizeof_val != sizeof(double))// cholmod only support double precision
            return 0;
        auto_ptr<cholmod_ctx> c(
            new cholmod_ctx(
                sizeof_int, sizeof_val, real_or_complex,
                nrows, ncols,
                ptr, idx, reinterpret_cast<const double *>(val), opts)
            );
        if(c.get() && !c->init(CHOLMOD_AUTO)){
            return c.release();
        }
    }
	// else if(!strcmp(name, "cholmod")) {
	// 	if(sizeof_val != sizeof(double))// cholmod only support double precision
	// 		return 0;
	// 	auto_ptr<cholmod_ctx> c(
	// 		new cholmod_ctx(
	// 			sizeof_int, sizeof_val, real_or_complex,
	// 			nrows, ncols,
	// 			ptr, idx, reinterpret_cast<const double *>(val), opts)
	// 		);
	// 	if(c.get() && !c->init())
	// 		return c.release();
	// }
	return 0;
}

int hj_sparse_direct_solver_A_solve(void *ctx, const void *b, void *x, size_t nrhs, void *opts)
{
	assert(ctx);
	solver_ctx *c = reinterpret_cast<solver_ctx *>(ctx);
	if(c->name_ == "umfpack") {
		umfpack_ctx *uc = reinterpret_cast<umfpack_ctx *>(ctx);;
	    return uc->solve(reinterpret_cast<const double *>(b), reinterpret_cast<double *>(x), nrhs, opts);
	}
	else if(c->name_ == "cholmod") {
		cholmod_ctx *c = reinterpret_cast<cholmod_ctx *>(ctx);;
	    return c->solve(reinterpret_cast<const double *>(b), reinterpret_cast<double *>(x), nrhs, opts);
	}
	return -1;
}

void hj_sparse_direct_solver_A_destroy(void *ctx)
{
	assert(ctx);
	solver_ctx *c = reinterpret_cast<solver_ctx *>(ctx);
	if(c->name_ == "umfpack") {
		umfpack_ctx *uc = reinterpret_cast<umfpack_ctx *>(ctx);
		delete uc;
	}
	else if(c->name_ == "cholmod") {
		cholmod_ctx *c = reinterpret_cast<cholmod_ctx *>(ctx);
		delete c;
	}
}

class laspack_ctx
{
public:
	laspack_ctx(size_t DIM)
		:A_(DIM, "A"), b_(DIM, "b"), x_(DIM, "x") {
	}
	laspack::QMatrix A_;
	laspack::Vector b_, x_;
};

void *hj_sparse_iterative_solver_A_create(
	unsigned char sizeof_int, 	unsigned char sizeof_val, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	const void *ptr, const void *idx, const void *val, //< link for internal use
	const char *name, //< copy to internal
	void *opts) //< link for internal use
{
	if(!name)
		name = "laspack";
	if(!strcmp(name, "laspack")) {
		if(nrows != ncols || sizeof_val != sizeof(double)) {
			cerr << "support quadric matrix in double only." << endl;
			return 0;
		}
		auto_ptr<laspack_ctx> c(new laspack_ctx(nrows));
		const double *pval = reinterpret_cast<const double *>(val);
		if(sizeof_int == sizeof(int)) {
			const int *pptr = reinterpret_cast<const int *>(ptr),
				*pidx = reinterpret_cast<const int *>(idx);
			c->A_.set(pptr, pidx, pval);
		}
		if(sizeof_int == sizeof(ptrdiff_t)) {
			const ptrdiff_t *pptr = reinterpret_cast<const ptrdiff_t *>(ptr),
				*pidx = reinterpret_cast<const ptrdiff_t *>(idx);
			c->A_.set(pptr, pidx, pval);
		}
		return c.release();
	}
	return 0;
}

int hj_sparse_iterative_solver_A_solve(void *ctx, const void *b, void *x, const size_t nrhs, void *opts)
{
	assert(ctx);
	laspack_ctx *c = reinterpret_cast<laspack_ctx *>(ctx);
	static const laspack_opts default_options("CG", "ILU");
	const laspack_opts *_opts = &default_options;
	if(opts)
		_opts = reinterpret_cast<laspack_opts *>(opts);
	IterProcType iter_proc = get_iter_proc(_opts->iter_name_);
	if(!iter_proc) {
		cerr << "no proper iterater solver found." << endl;
		return __LINE__;
	}
	PrecondProcType precond_proc = 0;
	if(_opts->precond_name_) {
		precond_proc = get_precond_proc(_opts->precond_name_);
		if(!precond_proc) {
			cerr << "no proper preconditioner found." << endl;
			return __LINE__;
		}
	}
	c->b_.set(reinterpret_cast<const double *>(b));
	c->x_.set(reinterpret_cast<const double *>(x));
	SetRTCAccuracy(_opts->precision_);
	iter_proc(&c->A_.get(), &c->x_.get(), &c->b_.get(), _opts->iter_num_, precond_proc, _opts->relax_);
	c->x_.get(reinterpret_cast<double *>(x));
	return 0;
}

void hj_sparse_iterative_solver_A_destroy(void *ctx)
{
	assert(ctx);
	laspack_ctx *c = reinterpret_cast<laspack_ctx *>(ctx);
	delete c;
}

