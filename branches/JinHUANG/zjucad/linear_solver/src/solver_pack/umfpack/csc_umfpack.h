#ifndef HJ_UMFPACK_WRAPPER_H_
#define HJ_UMFPACK_WRAPPER_H_

#include <stdlib.h>

#ifndef SUITESPARSE_WRAPPER_API
#define SUITESPARSE_WRAPPER_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

//! umfpack only support double precision

SUITESPARSE_WRAPPER_API
void *hj_umfpack_symbolic(
	unsigned char sizeof_int, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	const void *col_ptr, const void *row_idx, const double *values,
	void *rtn);

SUITESPARSE_WRAPPER_API
void *hj_umfpack_numeric(
	unsigned char sizeof_int, unsigned char real_or_complex,
	const void *sym,
	const void *col_ptr, const void *row_idx, const double *values,
	void *rtn);

//! return internal error code of UMFPACK
SUITESPARSE_WRAPPER_API
int hj_umfpack_solve(
	unsigned char sizeof_int, unsigned char real_or_complex,
	const void *num,
	const void *col_ptr, const void *row_idx, const double *values,
	const void *b, void *x);

SUITESPARSE_WRAPPER_API
void hj_umfpack_free_numeric(
	unsigned char sizeof_int, unsigned char real_or_complex,
	void **num);

SUITESPARSE_WRAPPER_API
void hj_umfpack_free_symbolic(
	unsigned char sizeof_int, unsigned char real_or_complex,
	void **sym);

// added by dzw

//step 0
void* csc_solver_new(unsigned char sizeof_int, unsigned char sizeof_val, unsigned char real_or_complex,
					 size_t nrows, size_t ncols,
					 const void *ptr, const void *idx, const void *val,
					 void *opts);

//step 1
int csc_reorder_matrix(void* solver_);

//step 2
int csc_analyze(void* solver_);

//step 3
int csc_factorize(void* solver_);

//step 4
int csc_solve(void* solver_, const void *b, void *x, size_t nrhs = 1, void *opts = 0);

//step 5
int csc_solver_delete(void* solver_);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

#include <cassert>
#include <memory>
#include <complex>

#include "hjlib/sparse/type_traits.h"

namespace hj { namespace sparse {

class nocopyable {
public:
	nocopyable(){}
private:
	nocopyable(const nocopyable&){}
	const nocopyable &operator =(const nocopyable&) { return *this; }
};

template <typename INT_TYPE = int, typename VAL_TYPE = double>
class umfpack
{
public:
	class symbolic : public nocopyable
	{
	public:
		~symbolic() {
			hj_umfpack_free_symbolic(sizeof(INT_TYPE), value_type<VAL_TYPE>::type,
									 &sym_);
		}

		symbolic(void *sym, const void *col_ptr, const void *row_idx)
			:sym_(sym), col_ptr_(col_ptr), row_idx_(row_idx) {
			assert(sym_ && col_ptr_ && row_idx_);
		}
		void *sym_;
		const void *col_ptr_, *row_idx_;
	};
	class numeric : public nocopyable
	{
	public:
		int solve(const VAL_TYPE *b, VAL_TYPE *x) const {
			return hj_umfpack_solve(sizeof(INT_TYPE), value_type<VAL_TYPE>::type,
									sym_.col_ptr_, sym_.row_idx_, values_,
									num_, b, x);
		}
		~numeric() {
			hj_umfpack_free_numeric(sizeof(INT_TYPE), value_type<VAL_TYPE>::type,
									&num_);
		}

		numeric(void *num, const symbolic &sym, const double *values)
			:num_(num), sym_(sym), values_(values) {
			assert(num_);
		}
		void *num_;
		const double *values_;
		const symbolic &sym_;
	};
	class factorization
	{
	public:
		int solve(const VAL_TYPE *b, VAL_TYPE *x) const {
			return num_->solve(b, x);
		}

		factorization(const symbolic *sym, const numeric *num)
			:sym_(sym), num_(num)
			{
				assert(sym && num);
			}

		const std::auto_ptr<const symbolic> sym_;
		const std::auto_ptr<const numeric> num_;
	};
};

template <typename INT_TYPE, typename VAL_TYPE>
typename umfpack<INT_TYPE, VAL_TYPE>::symbolic*
create_symbolic(
	INT_TYPE nrows, INT_TYPE ncols,
	const INT_TYPE *col_ptr, const INT_TYPE *row_idx, const VAL_TYPE *values = 0,
	INT_TYPE *rtn = 0) {
	void *sym = hj_umfpack_symbolic(sizeof(INT_TYPE), value_type<VAL_TYPE>::type,
									&nrows, &ncols, col_ptr, row_idx, values, rtn);
	if(sym)
		return new typename umfpack<INT_TYPE, VAL_TYPE>::symbolic(sym, col_ptr, row_idx);
	return 0;
}

template <typename INT_TYPE, typename VAL_TYPE>
typename umfpack<INT_TYPE, VAL_TYPE>::numeric*
create_numeric(
	const typename umfpack<INT_TYPE, VAL_TYPE>::symbolic &sym,
	const VAL_TYPE *values, INT_TYPE *rtn = 0) {
	void *num = hj_umfpack_numeric(sizeof(INT_TYPE), value_type<VAL_TYPE>::type,
								   sym.col_ptr_, sym.row_idx_,
								   sym.sym_, values, rtn);
	if(num)
		return new typename umfpack<INT_TYPE, VAL_TYPE>::numeric(num, sym, values);
	return 0;
}

template <typename INT_TYPE, typename VAL_TYPE>
typename umfpack<INT_TYPE, VAL_TYPE>::factorization *create_factorization(
	INT_TYPE nrows, INT_TYPE ncols,
	const INT_TYPE *col_ptr, const INT_TYPE *row_idx, const VAL_TYPE *values,
	INT_TYPE *rtn = 0) {
	std::auto_ptr<typename umfpack<INT_TYPE, VAL_TYPE>::symbolic> sym(
		create_symbolic(nrows, ncols, col_ptr, row_idx, values, rtn));
	if(!sym.get())
		return 0;
	std::auto_ptr<typename umfpack<INT_TYPE, VAL_TYPE>::numeric> num(
		create_numeric(*sym, values, rtn));
	if(!num.get())
		return 0;
	return new typename umfpack<INT_TYPE, VAL_TYPE>::factorization(sym.release(), num.release());
}

//! Typical usage 1:
//! auto_ptr<umfpack<int, double>::symbolic> sym = create_symbolic(...);
//! auto_ptr<umfpack<int, double>::numeric> num = create_numeric(...);
//! num->solve(...)

//! Typical usage 2:
//! auto_ptr<umfpack<>::factorization> fac = create_factorization(....)
//! num->solve(...)

}}

#endif

#endif
