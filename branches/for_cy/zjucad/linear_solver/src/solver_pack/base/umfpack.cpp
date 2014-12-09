extern "C" {
#include <suitesparse/umfpack.h>
}

#include <iostream>
using namespace std;

#include "./umfpack.h"

void *hj_umfpack_symbolic(
	unsigned char sizeof_int, unsigned char real_or_complex,
	size_t nrows, size_t ncols,
	const void *col_ptr, const void *row_idx, const double *values,
	void *rtn)
{
	void *sym = 0;
	if(sizeof_int == sizeof(int)) {
		const int nr = nrows, nc = ncols,
			*ptr = reinterpret_cast<const int *>(col_ptr),
			*idx = reinterpret_cast<const int *>(row_idx);
		int r = -1;
		if(real_or_complex == 'r')
			r = umfpack_di_symbolic(nr, nc, ptr, idx, values, &sym, 0, 0);
		else if(real_or_complex == 'c')
			r = umfpack_zi_symbolic(nr, nc, ptr, idx, values, 0, &sym, 0, 0);
		if(rtn)
			*reinterpret_cast<int *>(rtn) = r;
		if(r != UMFPACK_OK)
			sym = 0;
	}
	else if(sizeof_int == sizeof(UF_long)) {
		const UF_long nr = nrows, nc = ncols,
			*ptr = reinterpret_cast<const UF_long *>(col_ptr),
			*idx = reinterpret_cast<const UF_long *>(row_idx);
		UF_long r = -1;
		if(real_or_complex == 'r') {
			r = umfpack_dl_symbolic(nr, nc, ptr, idx, values, &sym, 0, 0);
		}
		else if(real_or_complex == 'c')
			r = umfpack_zl_symbolic(nr, nc, ptr, idx, values, 0, &sym, 0, 0);
		if(rtn)
			*reinterpret_cast<UF_long *>(rtn) = r;
		if(r != UMFPACK_OK)
			sym = 0;
	}

	return sym;
}

void *hj_umfpack_numeric(
	unsigned char sizeof_int, unsigned char real_or_complex,
	const void *sym,
	const void *col_ptr, const void *row_idx, const double *values,
	void *rtn)
{
	void *num = 0;
	if(sizeof_int == sizeof(int)) {
		const int *ptr = reinterpret_cast<const int *>(col_ptr),
			*idx = reinterpret_cast<const int *>(row_idx);
		int r = -1;
		if(real_or_complex == 'r')
			r = umfpack_di_numeric(ptr, idx, values, const_cast<void *>(sym), &num, 0, 0);
		else if(real_or_complex == 'c')
			r = umfpack_zi_numeric(ptr, idx, values, 0, const_cast<void*>(sym), &num, 0, 0);
		if(rtn)
			*reinterpret_cast<int *>(rtn) = r;
		if(r != UMFPACK_OK)
			num = 0;
	}
	else if(sizeof_int == sizeof(UF_long)) {
		const UF_long *ptr = reinterpret_cast<const UF_long *>(col_ptr),
			*idx = reinterpret_cast<const UF_long *>(row_idx);
		UF_long r = -1;
		if(real_or_complex == 'r')
			r = umfpack_dl_numeric(ptr, idx, values, const_cast<void*>(sym), &num, 0, 0);
		else if(real_or_complex == 'c')
			r = umfpack_zl_numeric(ptr, idx, values, 0, const_cast<void*>(sym), &num, 0, 0);
		if(rtn)
			*reinterpret_cast<UF_long *>(rtn) = r;
		if(r != UMFPACK_OK)
			num = 0;
	}

	return num;
}

int hj_umfpack_solve(
	unsigned char sizeof_int, unsigned char real_or_complex,
	const void *num,
	const void *col_ptr, const void *row_idx, const double *values,
	const void *b, void *x)
{
	const double *b0 = reinterpret_cast<const double *>(b);
	double *x0 = reinterpret_cast<double *>(x);
	if(sizeof_int == sizeof(int)) {
		const int *ptr = reinterpret_cast<const int *>(col_ptr),
			*idx = reinterpret_cast<const int *>(row_idx);
		if(real_or_complex == 'r')
			return umfpack_di_solve(UMFPACK_A, ptr, idx, values,
								 x0, b0, const_cast<void *>(num), 0, 0);
		else if(real_or_complex == 'c')
			return umfpack_zi_solve(UMFPACK_A, ptr, idx, values, 0,
								 x0, 0, b0, 0, const_cast<void *>(num), 0, 0);
	}
	else if(sizeof_int == sizeof(UF_long)) {
		const UF_long *ptr = reinterpret_cast<const UF_long *>(col_ptr),
			*idx = reinterpret_cast<const UF_long *>(row_idx);
		if(real_or_complex == 'r') {
			return umfpack_dl_solve(UMFPACK_A, ptr, idx, values,
								 x0, b0, const_cast<void *>(num), 0, 0);
		}
		else if(real_or_complex == 'c')
			return umfpack_zl_solve(UMFPACK_A, ptr, idx, values, 0,
								 x0, 0, b0, 0, const_cast<void *>(num), 0, 0);
	}
	return -1;
}

void hj_umfpack_free_numeric(
	unsigned char sizeof_int, unsigned char real_or_complex,
	void **num)
{
	if(sizeof_int == sizeof(int)) {
		if(real_or_complex == 'r')
			umfpack_di_free_numeric(num);
		else if(real_or_complex == 'c')
			umfpack_zi_free_numeric(num);
	}
	else if(sizeof_int == sizeof(UF_long)) {
		if(real_or_complex == 'r')
			umfpack_dl_free_numeric(num);
		else if(real_or_complex == 'c')
			umfpack_zl_free_numeric(num);
	}
}

void hj_umfpack_free_symbolic(
	unsigned char sizeof_int, unsigned char real_or_complex,
	void **sym)
{
	if(sizeof_int == sizeof(int)) {
		if(real_or_complex == 'r')
			umfpack_di_free_symbolic(sym);
		else if(real_or_complex == 'c')
			umfpack_zi_free_symbolic(sym);
	}
	else if(sizeof_int == sizeof(UF_long)) {
		if(real_or_complex == 'r')
			umfpack_dl_free_symbolic(sym);
		else if(real_or_complex == 'c')
			umfpack_zl_free_symbolic(sym);
	}
}
