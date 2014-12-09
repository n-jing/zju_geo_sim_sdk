#ifndef HJ_SPARSE_SOLVER_CTX_H_
#define HJ_SPARSE_SOLVER_CTX_H_

#include <inttypes.h>

#include <string>
#include <algorithm>
#include <iostream>
using namespace std;

#include "../umfpack/csc_umfpack.h"
#include "../cholmod/csc_cholmod.h"

extern "C" {
#include <suitesparse/cholmod.h>
}


class solver_ctx
{
public:
	solver_ctx(
		const char *name,
		unsigned char sizeof_int, unsigned char sizeof_val, unsigned char real_or_complex,
		size_t nrows, size_t ncols,
		const void *ptr, const void *idx, const double *val,
		void *opts)
		:name_(name),
		 sizeof_int_(sizeof_int), sizeof_val_(sizeof_val), real_or_complex_(real_or_complex),
		 nrows_(nrows), ncols_(ncols),
		 ptr_(ptr), idx_(idx), val_(val),
		 opts_(opts) {
	}
	virtual ~solver_ctx() {
	}
	std::string name_;
	unsigned char sizeof_int_, sizeof_val_, real_or_complex_;
	const size_t nrows_, ncols_;
	const void *ptr_, *idx_;
	const double *val_;
	void *opts_;
};

inline size_t get_size_t(unsigned char sizeof_int, const double *val)
{
	assert(val);
	if(sizeof_int == sizeof(uint32_t)) {
		return *reinterpret_cast<const uint32_t *>(val);
	}
	else if(sizeof_int == sizeof(uint64_t)) {
		return *reinterpret_cast<const uint64_t *>(val);
	}
	assert(0);
}

class umfpack_ctx : public solver_ctx
{
public:
	umfpack_ctx(
		unsigned char sizeof_int, unsigned char real_or_complex,
		size_t nrows, size_t ncols,
		const void *ptr, const void *idx, const double *val,
		void *opts)
		:solver_ctx("umfpack",
					sizeof_int, sizeof(double), real_or_complex,
					nrows, ncols, ptr, idx, val, opts),
		 sym_(0), num_(0) {
	}
	inline int init(void) {
		sym_ = hj_umfpack_symbolic(
			sizeof_int_, real_or_complex_,
			nrows_, ncols_, ptr_, idx_, val_, opts_);
		if(!sym_)
			return __LINE__;
		num_ = hj_umfpack_numeric(
			sizeof_int_, real_or_complex_,
			sym_,
			ptr_, idx_, val_,
			opts_);
		if(!num_)
			return __LINE__;
		return 0;
	}

	inline int solve(const double *b, double *x, size_t nrhs = 1, void *opts = 0) {
		assert(num_);
		const size_t stride = ((real_or_complex_=='r')?1:2)*nrows_;
		for(size_t i = 0; i < nrhs; ++i) {
			if(hj_umfpack_solve(
				sizeof_int_, real_or_complex_,
				num_,
				ptr_, idx_, val_,
				b, x))
				return __LINE__;
			b += stride;
			x += stride;
		}
		return 0;
	}

	~umfpack_ctx() {
		if(num_)
			hj_umfpack_free_numeric(
				sizeof_int_, real_or_complex_,
				&num_);
		if(sym_)
			hj_umfpack_free_symbolic(
				sizeof_int_, real_or_complex_,
				&sym_);
	}
private:
	void *sym_, *num_;
};

class cholmod_ctx : public solver_ctx
{
public:
	cholmod_ctx(
		unsigned char sizeof_int, unsigned char sizeof_val, unsigned char real_or_complex,
		size_t nrows, size_t ncols,
		const void *ptr, const void *idx, const double *val,
		void *opts)
		:solver_ctx("cholmod",
					sizeof_int, sizeof_val, real_or_complex,
					nrows, ncols, ptr, idx, val, opts),
		 L_(0)
		{
			c_.ptr = 0;
			A_.ptr = 0;
			A_ = hj_link_cholmod_sparse_to_csc(
				sizeof_int_, real_or_complex_,
				nrows_, ncols_,
				const_cast<void *>(ptr_), const_cast<void *>(idx_), const_cast<double *>(val_),
				1, 1);
	}
      inline int init(int method = CHOLMOD_SUPERNODAL) {
		if(!hj_cholmod_start(sizeof_int_, &c_))
			return __LINE__;
        ((cholmod_common*)c_.ptr)->supernodal = method;
		L_ = hj_cholmod_analyze(&A_, &c_);
		if(!L_->ptr)
			return __LINE__;
		if(hj_cholmod_factorize(&A_, L_, &c_))
			return __LINE__;
		return 0;
	}

	/* inline int init(void) { */
	/* 	if(!hj_cholmod_start(sizeof_int_, &c_)) */
	/* 		return __LINE__; */
	/* 	L_ = hj_cholmod_analyze(&A_, &c_); */
	/* 	if(!L_->ptr) */
	/* 		return __LINE__; */
	/* 	if(hj_cholmod_factorize(&A_, L_, &c_)) */
	/* 		return __LINE__; */
	/* 	return 0; */
	/* } */

	inline int solve(const double *b, double *x, size_t nrhs = 1, void *opts = 0) {
		hj_cholmod_dense B = hj_link_cholmod_dense(
			real_or_complex_,
			nrows_, const_cast<double *>(b), nrhs, nrows_);
		hj_cholmod_dense *X = hj_cholmod_solve(0, L_, &B, &c_);
		if(!X->ptr)
			return __LINE__;
		hj_copy_out_cholmod_dense(X, x, 0);
		hj_cholmod_free_dense(&X, &c_);
		hj_unlink_cholmod_dense(&B);
		return 0;
	}

	~cholmod_ctx() {
		if(L_)
			hj_cholmod_free_factor(&L_, &c_);
		if(c_.ptr)
			hj_cholmod_finish(&c_);
		if(A_.ptr)
			hj_unlink_cholmod_sparse_to_csc(&A_);
	}
private:
	hj_cholmod_common c_;
	hj_cholmod_sparse A_;
	hj_cholmod_factor *L_;
};

#endif
