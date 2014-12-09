extern "C" {
#include <suitesparse/umfpack.h>
}

#include <iostream>
using namespace std;

#include "./csc_umfpack.h"

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



/////////////////////////////////////added by dzw
class umfpack_ctx
{
public:
	umfpack_ctx(
		unsigned char sizeof_int, unsigned char sizeof_val, unsigned char real_or_complex,
		size_t nrows, size_t ncols,
		const void *ptr, const void *idx, const double *val,
		void *opts)
	  :sizeof_int_(sizeof_int), sizeof_val_(sizeof(double)), real_or_complex_(real_or_complex),
	        nrows_(nrows), ncols_(ncols),
		ptr_(ptr), idx_(idx), val_(val),
		opts_(opts),
		sym_(0), num_(0) {
	}
	int umfpack_symbolic() {
		sym_ = hj_umfpack_symbolic(
			sizeof_int_, real_or_complex_,
			nrows_, ncols_, ptr_, idx_, val_, opts_);
		if(!sym_)
			return __LINE__;
		return 0;
	}
	int umfpack_numeric(){
		num_ = hj_umfpack_numeric(
			sizeof_int_, real_or_complex_,
			sym_,
			ptr_, idx_, val_,
			opts_);
		if(!num_)
			return __LINE__;
		return 0;
	}

	int solve(const double *b, double *x, size_t nrhs = 1, void *opts = 0) {
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
	if(sizeof_val != sizeof(double))// umfpack only support double precision
		return 0;
	auto_ptr<umfpack_ctx> umfpack_solver(new umfpack_ctx(sizeof_int,sizeof_val,real_or_complex,nrows,ncols,ptr,idx,reinterpret_cast<const double *>(val),opts));
	return umfpack_solver.release();
}

int csc_reorder_matrix(void* solver_){
	return 0;
}

//step 2
int csc_analyze(void* solver_){
	umfpack_ctx *umfpack_solver = reinterpret_cast<umfpack_ctx*>(solver_);
	return umfpack_solver->umfpack_symbolic();
}

//step 3
int csc_factorize(void* solver_){
	umfpack_ctx *umfpack_solver = reinterpret_cast<umfpack_ctx*>(solver_);
	return umfpack_solver->umfpack_numeric();
}

//step 4
int csc_solve(void* solver_, const void *b, void *x, size_t nrhs, void *opts){
	umfpack_ctx *umfpack_solver = reinterpret_cast<umfpack_ctx*>(solver_);
	return umfpack_solver->solve(reinterpret_cast<const double *>(b),reinterpret_cast<double *>(x),nrhs,opts);
}

//step 5
int csc_solver_delete(void* solver_){
	assert(solver_);
	umfpack_ctx *umfpack_solver = reinterpret_cast<umfpack_ctx*>(solver_);
	delete umfpack_solver; 
	return 0;
}
