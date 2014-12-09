#include "./cholmod.h"

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
		cholmod_factorize(pA, pL, pc);
	else if(pc->itype == CHOLMOD_LONG)
		cholmod_l_factorize(pA, pL, pc);
	else
		return -1;
  return pc->status == CHOLMOD_NOT_POSDEF;
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
