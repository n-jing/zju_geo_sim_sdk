#ifndef HJ_FAST_AAT_H_
#define HJ_FAST_AAT_H_

#include <stdint.h>

#include <algorithm>
#include <iostream>

#include <hjlib/sparse/sparse.h>

//! @param AAT must be a sorted csc with correct patten

//! @param method: which_part: -1 lower, 0 full, 1 upper, -1 is not
//! supported currently.

template <typename CSC>
void fast_AAT(const CSC &A,
			  CSC &AAT,
			  bool is_A_sorted = false,
			  int which_part = 0)
{
  typedef typename CSC::val_type val_type;
  typedef typename CSC::int_type int_type;

  using namespace std;
  using namespace hj::sparse;
  using namespace zjucad::matrix;

	if(which_part == -1) {
		cerr << "not support currently." << endl;
		return;
	}
	fill(&AAT.val()[0], &AAT.val()[0]+AAT.nnz(), 0);
	assert(A.size(1) == AAT.size(1) && A.size(1) == AAT.size(2));
	const int_type n = A.size(2);
	int_type i, col_i = 0, row_i = 0;
	for(i = 0; i < n; ++i) {
		for(col_i = A.ptr()[i]; col_i < A.ptr()[i+1]; ++col_i) {	// for each nz column
			const int_type col_idx_of_AAT = A.idx()[col_i];
			const int_type nz_beg = AAT.ptr()[col_idx_of_AAT], nz_end = AAT.ptr()[col_idx_of_AAT+1];
			for(row_i = A.ptr()[i]; row_i < A.ptr()[i+1]; ++row_i) {	// scalar * sparse, for each nz row
				if(A.idx()[row_i] > col_idx_of_AAT) {
					if(is_A_sorted) break;
					else continue;
				}
				const int_type nz_of_AAT =
					lower_bound(&AAT.idx()[nz_beg], &AAT.idx()[0]+nz_end, A.idx()[row_i])
					-&AAT.idx()[0];
				AAT.val()[nz_of_AAT] += A.val()[col_i]*A.val()[row_i];
			}
		}
	}

	// get another part
	for(col_i = 0; col_i < AAT.size(2) && which_part == 0; ++col_i) {
		for(size_t nzi = AAT.ptr()[col_i]; nzi < AAT.ptr()[col_i+1]; ++nzi) {
			const size_t row_i = AAT.idx()[nzi];
			if(row_i > col_i) {
				if(AAT.val()[nzi] != 0) {
					cerr << "error: " << AAT.val()[nzi] << endl;
				}
				const size_t offset = lower_bound(&AAT.idx()[AAT.ptr()[row_i]], &AAT.idx()[0]+AAT.ptr()[row_i+1], col_i)-&AAT.idx()[0];
				AAT.val()[nzi] = AAT.val()[offset];
			}
		}
	}
}

#endif
