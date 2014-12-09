#include "./csc_cholmod.h"

extern "C" {
#include <suitesparse/cholmod.h>
}

#include <cassert>
#include <memory>
#include <iostream>

#include <string.h>

using namespace std;

#define CHOLMOD(type, val) \
	reinterpret_cast<cholmod_##type *>(val->ptr)
#define REF(type, val) \
	cholmod_##type *p##val = reinterpret_cast<cholmod_##type *>(val->ptr);

int hj_cholmod_start(unsigned char sizeof_int, hj_cholmod_common *c)
{
	cholmod_common *cc = new cholmod_common;
	int rtn = -1;
	if(sizeof_int == sizeof(int))
		rtn = cholmod_start(cc);
	else if(sizeof_int == sizeof(UF_long))
		rtn = cholmod_l_start(cc);
	c->ptr = cc;
	return rtn;
}

int hj_cholmod_finish(hj_cholmod_common *c)
{
	REF(common, c);
	int rtn = -1;
	if(pc->itype == CHOLMOD_INT)
		rtn = cholmod_finish(pc);
	else if(pc->itype == CHOLMOD_LONG)
		rtn = cholmod_l_finish(pc);
	delete pc;
	return rtn;
}

// hj_cholmod_sparse hj_cholmod_aat(hj_cholmod_sparse *A, void *fset, void *fsize, int mode, hj_cholmod_common *c)
// {
// 	REF(common, c);
// 	REF(sparse, A);
// 	hj_cholmod_sparse rtn;
// 	rtn. ptr = 0;
// 	if(pc->itype == CHOLMOD_INT)
// 		rtn.ptr = cholmod_aat(pA, reinterpret_cast<int *>(fset), *reinterpret_cast<int *>(fsize), mode, pc);
// 	else if(pc->itype == CHOLMOD_LONG)
// 		rtn.ptr = cholmod_l_aat(CHOLMOD(sparse, A), reinterpret_cast<UF_long *>(fset), *reinterpret_cast<UF_long *>(fsize), mode, pc);
// 	else
// 		assert(0);
// 	return rtn;
// }

hj_cholmod_factor *hj_cholmod_analyze(hj_cholmod_sparse *A, hj_cholmod_common *c)
{
	REF(common, c);
	REF(sparse, A);
	auto_ptr<hj_cholmod_factor> rtn(new hj_cholmod_factor);
	rtn->ptr = 0;
	if(pc->itype == CHOLMOD_INT)
		rtn->ptr = cholmod_analyze(pA, pc);
	else if(pc->itype == CHOLMOD_LONG)
		rtn->ptr = cholmod_l_analyze(pA, pc);
	else
		assert(0);
	return rtn.release();
}

int hj_cholmod_factorize(hj_cholmod_sparse *A, hj_cholmod_factor *L, hj_cholmod_common *c)
{
	REF(common, c);
	REF(sparse, A);
	REF(factor, L);
	if(pc->itype == CHOLMOD_INT)
		return cholmod_factorize(pA, pL, pc);
	else if(pc->itype == CHOLMOD_LONG)
		return cholmod_l_factorize(pA, pL, pc);
	else
		return -1;
}

int hj_cholmod_free_factor(hj_cholmod_factor **L, hj_cholmod_common *c)
{
	REF(common, c);
	cholmod_factor **ppL = reinterpret_cast<cholmod_factor **>(&(*L)->ptr);
	int rtn = -1;
	if(pc->itype == CHOLMOD_INT)
		rtn = cholmod_free_factor(ppL, pc);
	else if(pc->itype == CHOLMOD_LONG)
		rtn = cholmod_l_free_factor(ppL, pc);
	delete *L;
	return rtn;
}

hj_cholmod_dense *hj_cholmod_solve(int sys, hj_cholmod_factor *L, hj_cholmod_dense *B, hj_cholmod_common *c)
{
	REF(common, c);
	REF(factor, L);
	REF(dense, B);
	auto_ptr<hj_cholmod_dense> rtn(new hj_cholmod_dense);
	rtn->ptr = 0;
	if(pc->itype == CHOLMOD_INT)
		rtn->ptr = cholmod_solve(sys, pL, pB, pc);
	else if(pc->itype == CHOLMOD_LONG)
		rtn->ptr = cholmod_l_solve(sys, pL, pB, pc);
	else
		assert(0);
	return rtn.release();
}

int hj_cholmod_free_dense(hj_cholmod_dense **X, hj_cholmod_common *c)
{
	REF(common, c);
	cholmod_dense **ppX = reinterpret_cast<cholmod_dense **>(&(*X)->ptr);
	int rtn = -1;
	if(pc->itype == CHOLMOD_INT)
		rtn = cholmod_free_dense(ppX, pc);
	else if(pc->itype == CHOLMOD_LONG)
		rtn = cholmod_l_free_dense(ppX, pc);
	delete *X;
	return rtn;
}

hj_cholmod_sparse hj_link_cholmod_sparse_to_csc(
	unsigned char sizeof_int, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	void *col_ptr, void *row_idx, double *values,
	int stype, int is_sorted)
{
	auto_ptr<cholmod_sparse> to(new cholmod_sparse);

	hj_cholmod_sparse rtn;
	rtn.ptr = 0;

	memset(to.get(), 0, sizeof(cholmod_sparse));

	to->nrow = nrows;
	to->ncol = ncols;
	if(sizeof_int == sizeof(int)) {
		const int *ptr = reinterpret_cast<const int *>(col_ptr);
		to->nzmax = ptr[to->ncol]-ptr[0];
	}
	else if(sizeof_int == sizeof(UF_long)) {
		const UF_long *ptr = reinterpret_cast<const UF_long *>(col_ptr);
		to->nzmax = ptr[to->ncol]-ptr[0];
	}
	else {
		return rtn;
	}

	to->p = col_ptr;
	to->i = row_idx;
	to->x = values;

	to->stype = stype;

	// indeed, cholmod support mixed integer type by CHOLMOD_INTLONG, which separate p from i and nz.
	to->itype = (sizeof_int==sizeof(int)?CHOLMOD_INT:CHOLMOD_LONG);
	to->xtype = ((real_or_complex=='r')?CHOLMOD_REAL:CHOLMOD_COMPLEX);
	to->dtype = CHOLMOD_DOUBLE;

	to->sorted = (is_sorted?1:0);
	to->packed = 1;

	rtn.ptr = to.release();
	return rtn;
}

void hj_unlink_cholmod_sparse_to_csc(hj_cholmod_sparse *A)
{
	REF(sparse, A);
	delete pA;
}

hj_cholmod_dense hj_link_cholmod_dense(
	unsigned char real_or_complex,
	size_t nrows, double *values,
	size_t ncols, size_t ld)
{
	auto_ptr<cholmod_dense> to(new cholmod_dense);
	hj_cholmod_dense rtn;
	rtn.ptr = 0;

	to->nrow = nrows;
	to->ncol = ncols;
	to->d = ld;
	to->nzmax = to->nrow*to->ncol;

	to->x = values;
	to->xtype = ((real_or_complex=='r')?CHOLMOD_REAL:CHOLMOD_COMPLEX);
	to->dtype = CHOLMOD_DOUBLE;

	rtn.ptr = to.release();
	return rtn;
}

void hj_unlink_cholmod_dense(hj_cholmod_dense *X)
{
	REF(dense, X);
	delete pX;
}

int hj_copy_out_cholmod_dense(
	const hj_cholmod_dense *X, double *values,
	const size_t *ld)
{
	const cholmod_dense *pX = reinterpret_cast<const cholmod_dense *>(X->ptr);

	if(pX->xtype == CHOLMOD_ZOMPLEX)
		return -1; // not support currently

	const size_t ld0 = (ld?*ld:pX->nrow),
		byte_per_value = sizeof(double)*((pX->xtype==CHOLMOD_COMPLEX)?2:1),
		chunk_size = byte_per_value*pX->nrow,
		stride[2] = {
		byte_per_value*ld0, byte_per_value*pX->d
	};
	for(size_t ci = 0; ci < pX->ncol; ++ci) {
		memcpy(reinterpret_cast<char *>(values)+ci*stride[0],
			   reinterpret_cast<const char *>(pX->x)+ci*stride[1],
			   chunk_size);
	}
	return 0;
}

//////////////////////////////////////////////////////added by dzw

class cholmod_ctx
{
public:
	cholmod_ctx(
		unsigned char sizeof_int, unsigned char sizeof_val, unsigned char real_or_complex,
		size_t nrows, size_t ncols,
		const void *ptr, const void *idx, const double *val,
		void *opts)
		:sizeof_int_(sizeof_int), sizeof_val_(sizeof_val), real_or_complex_(real_or_complex),
		nrows_(nrows), ncols_(ncols),
		ptr_(ptr), idx_(idx), val_(val),
		opts_(opts),
		L_(0)
	{
		c_.ptr = 0;
		A_.ptr = 0;
                // convert csc(ptr, idx, val) to cholmod_sparse
		A_ = hj_link_cholmod_sparse_to_csc(
			sizeof_int_, real_or_complex_,
			nrows_, ncols_,
			const_cast<void *>(ptr_), const_cast<void *>(idx_), const_cast<double *>(val_),
			1, 1);
	}
	int cholmod_start(){
		if(!hj_cholmod_start(sizeof_int_, &c_))
			return __LINE__;
		return 0;
	}
	int cholmod_analyze(){
		L_ = hj_cholmod_analyze(&A_, &c_);
		if(!L_->ptr)
			return __LINE__;
		return 0;
	}
	int cholmod_factorize(){
		if(!hj_cholmod_factorize(&A_, L_, &c_))
			return __LINE__;
		return 0;
	}
	int solve(const double *b, double *x, size_t nrhs = 1, void *opts = 0) {
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
	unsigned char sizeof_int_, sizeof_val_, real_or_complex_;
	const size_t nrows_, ncols_;
	const void *ptr_, *idx_;
	const double *val_;
	void *opts_;
};

void* csc_solver_new(unsigned char sizeof_int, unsigned char sizeof_val, unsigned char real_or_complex,
					size_t nrows, size_t ncols,
					const void *ptr, const void *idx, const void *val,
					void *opts)
{
	if(sizeof_val != sizeof(double))// cholmod only support double precision
		return 0;
	auto_ptr<cholmod_ctx> cholmod_solver(new cholmod_ctx(sizeof_int,sizeof_val,real_or_complex,nrows,ncols,ptr,idx,reinterpret_cast<const double *>(val),opts));
	return cholmod_solver.release();
}

int csc_reorder_matrix(void* solver_){
	cholmod_ctx *cholmod_solver = reinterpret_cast<cholmod_ctx*>(solver_);
	return cholmod_solver->cholmod_start();
}

//step 2
int csc_analyze(void* solver_){
	cholmod_ctx *cholmod_solver = reinterpret_cast<cholmod_ctx*>(solver_);
	return cholmod_solver->cholmod_analyze();
}

//step 3
int csc_factorize(void* solver_){
	cholmod_ctx *cholmod_solver = reinterpret_cast<cholmod_ctx*>(solver_);
	return cholmod_solver->cholmod_factorize();
}

//step 4
int csc_solve(void* solver_, const void *b, void *x, size_t nrhs, void *opts){
	cholmod_ctx *cholmod_solver = reinterpret_cast<cholmod_ctx*>(solver_);
	return cholmod_solver->solve(reinterpret_cast<const double *>(b),reinterpret_cast<double *>(x),nrhs,opts);
}

//step 5
int csc_solver_delete(void* solver_){
	cholmod_ctx *cholmod_solver = reinterpret_cast<cholmod_ctx*>(solver_);
	delete cholmod_solver; 
	return 0;
}
