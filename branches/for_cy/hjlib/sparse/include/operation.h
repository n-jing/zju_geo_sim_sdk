#ifndef HJ_SPARSE_OPERATION_H_
#define HJ_SPARSE_OPERATION_H_

#include "format.h"

namespace hj { namespace sparse {

//! @brief B += A
//! NOTICE: B must be pre-allocated!!
template<
	typename T1, typename INT_TYPE1,
	typename T2, typename INT_TYPE2,
	template <typename K, typename T> class MAP_TYPE>
csc_by_vm<T2, INT_TYPE2, MAP_TYPE> &
acc(const csc<T1, INT_TYPE1> &A, csc_by_vm<T2, INT_TYPE2, MAP_TYPE> &B)
{
	INT_TYPE1 ci;
	for(ci = 0; ci < A.size(2); ++ci) {
		for(INT_TYPE1 ri = A.ptr()[ci]; ri < A.ptr()[ci+1]; ++ri) {
			B[ci][A.idx()[ri]] += A.val()[ri];
		}
	}
	return B;
}

//! NOTICE: B's pattern must equal or larger than A
template<
	typename T1, typename INT_TYPE1,
	typename T2, typename INT_TYPE2>
csc<T2, INT_TYPE2> &
acc(const csc<T1, INT_TYPE1> &A, csc<T2, INT_TYPE2> &B)
{
	INT_TYPE1 ci;
	for(ci = 0; ci < A.size(2); ++ci) { // for each col
		for(INT_TYPE1 ri = A.ptr()[ci]; ri < A.ptr()[ci+1]; ++ri) {
			for(INT_TYPE2 Bri = B.ptr()[ci]; Bri < B.ptr()[ci+1]; ++Bri) {
				if(A.idx()[ri] == B.idx()[Bri]) {
					B.val()[Bri] += A.val()[ri];
					break;
				}
			}
		}
	}
	return B;
}

//! NOTICE: B's pattern must equal or larger than A and is SORTED!
template<
	typename T1, typename INT_TYPE1,
	typename T2, typename INT_TYPE2>
csc<T2, INT_TYPE2> &
acc_sorted(const csc<T1, INT_TYPE1> &A, csc<T2, INT_TYPE2> &B)
{
	INT_TYPE1 ci;
	for(ci = 0; ci < A.size(2); ++ci) { // for each col
		for(INT_TYPE1 ri = A.ptr()[ci]; ri < A.ptr()[ci+1]; ++ri) {
			const INT_TYPE2 *pos = std::lower_bound(
				&B.idx()[0]+B.ptr()[ci],
				&B.idx()[0]+B.ptr()[ci+1],
				A.idx()[ri]);
			B.val()[pos - &B.idx()[0]] += A.val()[ri];
		}
	}
	return B;
}

//! @brief NOTICE: C += A*B
template <typename T, typename INT_TYPE, typename M1, typename M2>
M2 &mm(bool transA, const csc<T, INT_TYPE> &A,
	   const zjucad::matrix::matrix_expression<M1> &B,
	   zjucad::matrix::matrix_expression<M2> &C)
{
	using namespace zjucad::matrix;
	assert(B().size(1) == (transA?A.size(1):A.size(2)));
	assert(C().size(1) == (transA?A.size(2):A.size(1))
		   && C().size(2) == B().size(2));
#define mm_OP(op) \
	for(INT_TYPE ci = 0; ci < A.size(2); ++ci) { for(; vi < A.ptr()[ci+1]; ++vi) { \
		op;\
	}}

	INT_TYPE vi = A.ptr()[0];
	if(transA) {
		mm_OP(C()(ci, colon()) += A.val()[vi]*B()(A.idx()[vi], colon()));
	}
	else {
		mm_OP(C()(A.idx()[vi], colon()) += A.val()[vi]*B()(ci, colon()));
	}
#undef mm_OP
	return C();
}

//! @brief y += A*x
template <typename CSC, typename V1, typename V2>
V2 &mv(bool transA, const CSC &A, const V1 &x, V2 &y)
{
#define mv_OP(op) \
	for(size_type ci = 0; ci < A.size(2); ++ci) { \
		for(; vi < A.ptr()[ci+1]; ++vi) {		  \
			op;									  \
		}}

	typename CSC::int_type vi = A.ptr()[0];
	if(transA) {
		mv_OP(y[ci] += A.val()[vi]*x[A.idx()[vi]]);
	}
	else {
		mv_OP(y[A.idx()[vi]] += A.val()[vi]*x[ci]);
	}
#undef mv_OP
	return y;
}

//! @brief transpose
template <typename T1, typename INT_TYPE1, typename T2, typename INT_TYPE2>
csc<T2, INT_TYPE2> &trans(const csc<T1, INT_TYPE1> &A,csc<T2, INT_TYPE2> &AT)
{
	INT_TYPE1 nz = nnz(A);
	AT.resize(A.size(2), A.size(1), nz);
	std::vector<INT_TYPE1> rowvec;
	rowvec.resize(A.size(1),0);
	for(INT_TYPE1 i=0;i<nz;++i)
		++rowvec[A.idx()[i]];
	AT.ptr()[0] = 0;
	AT.ptr()[1] = 0;
	for(INT_TYPE1 i=1;i<A.size(1);++i)
		AT.ptr()[i+1] = AT.ptr()[i]+rowvec[i-1];
	INT_TYPE1 pos;
	for(INT_TYPE1 i = 0; i < A.size(2); ++i) {  
		for(pos = A.ptr()[i]; pos < A.ptr()[i+1]; ++pos) { 
			INT_TYPE1 ind = AT.ptr()[A.idx()[pos]+1]++;
			AT.idx()[ind] = i;
			AT.val()[ind] = A.val()[pos];
		}
	}
	assert(AT.ptr()[A.size(1)]==nz);
	return AT;
}

//! @brief compute AAT += A*A'
//! NOTICE AAT must be preallocated
template <typename T1, typename INT_TYPE1,
		  typename T2, typename INT_TYPE2,
		  template <typename K, typename T> class MAP_TYPE>
csc_by_vm<T2, INT_TYPE2, MAP_TYPE> &
acc_AAT(const csc<T1, INT_TYPE1> &A, csc_by_vm<T2, INT_TYPE2, MAP_TYPE> &AAT)
{
	assert(A.size(1) == AAT.size(1) && A.size(1) == AAT.size(2));
	const INT_TYPE2 n = A.size(2);
	INT_TYPE2 i, col_i = 0, row_i = 0;
	for(i = 0; i < n; ++i) {	// vector outer product: vvT(i, j) = vi*vj, for each v
		for(col_i = A.ptr()[i]; col_i < A.ptr()[i+1]; ++col_i) {	// for each nz column
      //			if(A.val()[col_i] == 0) continue;
			map_vec<T2, INT_TYPE2, MAP_TYPE> &col = AAT[A.idx()[col_i]];
			for(row_i = A.ptr()[i]; row_i < A.ptr()[i+1]; ++row_i)	// scalar * sparse, for each nz row
				col[A.idx()[row_i]] += A.val()[col_i]*A.val()[row_i];
		}
	}
	return AAT;
}

//! @brief compute AAT = A*A'
template<template <typename K, typename T> class MAP_TYPE = map_by_unsorted_vector>
class AAT
{
public:
	template<typename T1, typename INT_TYPE1, typename T2, typename INT_TYPE2>
	AAT(const csc<T1, INT_TYPE1> &A, csc<T2, INT_TYPE2> &AAT) {
		csc_by_vm<T2, INT_TYPE2, MAP_TYPE> vm(A.size(1), A.size(1), 0); // NOTICE nz per col
		acc_AAT(A, vm);
		convert(vm, AAT);
	}
};

// template <template <typename K, typename T> class MAP_TYPE,
// 		  typename T1, typename INT_TYPE1, typename T2, typename INT_TYPE2>
// csc<T2, INT_TYPE2> &AAT(const csc<T1, INT_TYPE1> &A, csc<T2, INT_TYPE2> &AAT)
// {
// 	std::vector<map_vec<T2, INT_TYPE2, MAP_TYPE> > vm(A.size(1));
// 	acc_AAT(A, vm);
// 	return convert(vm, AAT, A.size(1));
// }


//! @brief computer C+=A*B'
//! NOTICE: C must be preallocated!
template <typename T1, typename INT_TYPE1,
		  typename T2, typename INT_TYPE2,
		  typename T3, typename INT_TYPE3,
		  template <typename K, typename T> class MAP_TYPE>
csc_by_vm<INT_TYPE3, T3, MAP_TYPE> &
acc_ABT(const csc<T1, INT_TYPE1> &A, const csc<T2, INT_TYPE2> &B, csc_by_vm<INT_TYPE3, T3, MAP_TYPE> &C)
{
	const INT_TYPE1 k = A.size(2);
	assert(k == B.size(2));
	assert(A.size(1) == C.size(1));
	assert(B.size(1) == C.size(2));

	INT_TYPE1 i, row_i = 0;
	INT_TYPE2 col_i = 0;
	for(i = 0; i < k; ++i) {    // vector outer product: vvT(i, j) = vi*vj, for each v
			for(col_i = B.ptr()[i]; col_i < B.ptr()[i+1]; ++col_i) {    // for each nz column
        //				if(B.val()[col_i] == 0) continue;
				map_vec<INT_TYPE3, T3, MAP_TYPE> &col = C[B.idx()[col_i]];
				for(row_i = A.ptr()[i]; row_i < A.ptr()[i+1]; ++row_i)    // scalar * sparse, for each nz row
					col[A.idx()[row_i]] += B.val()[col_i]*A.val()[row_i];
			}
	}
	return C;
}

template<template <typename K, typename T> class MAP_TYPE = map_by_unsorted_vector>
class ABT {
public:
	template <typename T1, typename INT_TYPE1,
			  typename T2, typename INT_TYPE2,
			  typename T3, typename INT_TYPE3>
	ABT(const csc<T1, INT_TYPE1> &A, const csc<T2, INT_TYPE2> &B, csc<T3, INT_TYPE3> &C) {
		csc_by_vm<T3, INT_TYPE3, MAP_TYPE> vm(A.size(1), B.size(1), 0);  // NOTICE nz per col
		acc_ABT(A, B, vm);
		convert(vm, C);
	}
};

template<template <typename K, typename T> class MAP_TYPE = map_by_unsorted_vector>
class MM {
public:
	//! the most efficient is ABT
	template <typename T1, typename INT_TYPE1,
			  typename T2, typename INT_TYPE2,
			  typename T3, typename INT_TYPE3>
	MM(bool transA, const csc<T1, INT_TYPE1> &A,
	   bool transB, const csc<T2, INT_TYPE2> &B,
	   csc<T3, INT_TYPE3> &C) {
		if(transA == true) {
			csc<T1, INT_TYPE1> AT;
			MM(false, trans(A, AT), transB, B, C);
			return;
		}
		if(transB == false) {
			csc<T1, INT_TYPE1> BT;
			MM(transA, A, true, trans(B, BT), C);
			return;
		}
		assert(transA == false && transB == true);
		ABT<MAP_TYPE>(A, B, C);
	}
};

template <typename VAL_TYPE, typename INT_TYPE>
bool is_sorted_csc(const csc<VAL_TYPE, INT_TYPE> &A)
{
	for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
		for(INT_TYPE nzi = A.ptr()[ci]+1; nzi < A.ptr()[ci+1]; ++nzi) {
			if(A.idx()[nzi-1] >= A.idx()[nzi])
				return false;
		}
	}
	return true;
}

}}
#endif
