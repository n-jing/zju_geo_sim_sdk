#ifndef HJ_SPARSE_CHOLMOD_H_
#define HJ_SPARSE_CHOLMOD_H_

#include <stdlib.h>

extern "C" {

//! cholmod only support double precision, see the manual

#define HJ_CHOLMOD_TYPE(name)  \
	struct hj_cholmod_##name { \
		void *ptr;			  \
	};

HJ_CHOLMOD_TYPE(common)
HJ_CHOLMOD_TYPE(sparse)
HJ_CHOLMOD_TYPE(factor)
HJ_CHOLMOD_TYPE(dense)

int hj_cholmod_start(unsigned char sizeof_int, hj_cholmod_common *c);
int hj_cholmod_finish(hj_cholmod_common *c);


// hj_cholmod_sparse hj_cholmod_aat(hj_cholmod_sparse *A, void *fset, void *fsize, int mode, hj_cholmod_common *c);

hj_cholmod_factor *hj_cholmod_analyze(hj_cholmod_sparse *A, hj_cholmod_common *c);
int hj_cholmod_factorize(hj_cholmod_sparse *A, hj_cholmod_factor *L, hj_cholmod_common *c);
int hj_cholmod_free_factor(hj_cholmod_factor **L, hj_cholmod_common *c);

//! @param sys 0:Ax=b
hj_cholmod_dense *hj_cholmod_solve(int sys, hj_cholmod_factor *L, hj_cholmod_dense *B, hj_cholmod_common *c);

int hj_cholmod_free_dense(hj_cholmod_dense **X, hj_cholmod_common *c);

//! convert

//! @param stype 0: "unsymmetric", >0: upper of symmetric, <0: lower of symmetric
hj_cholmod_sparse hj_link_cholmod_sparse_to_csc(
	unsigned char sizeof_int, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	void *col_ptr, void *row_idx, double *values,
	int stype, int is_sorted);

void hj_unlink_cholmod_sparse_to_csc(hj_cholmod_sparse *A);

hj_cholmod_dense hj_link_cholmod_dense(
	unsigned char real_or_complex,
	size_t nrows, double *values,
	size_t ncols, size_t ld
);
	
void hj_unlink_cholmod_dense(hj_cholmod_dense *X);

int hj_copy_out_cholmod_dense(
	const hj_cholmod_dense *A, double *values,
	const size_t *ld //< optional
	);

}


#undef HJ_CHOLMOD_TYPE
#endif
