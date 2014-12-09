#ifndef HJ_SPARSE_LASPACK_H_
#define HJ_SPARSE_LASPACK_H_

#include <cassert>

extern "C" {
#include <laspack/qmatrix.h>
#include <laspack/rtc.h>
#include <laspack/itersolv.h>
#include <laspack/precond.h>
}

#include "include/sparse.h"

namespace laspack {

//! NOTICE: laspack index of matrix and vector start from 1 instead of
//! 0, but the entry offset start from 0.

class QMatrix
{
public:
	QMatrix(size_t DIM, const char *name) {
		Q_Constr(&M_, const_cast<char *>(name), DIM, False, Clmws, Normal, True);
	}
	template <typename VAL_TYPE, typename INT_TYPE>
	void set(const INT_TYPE *ptr, const INT_TYPE *idx, const VAL_TYPE *val) {
		//! TODO: strange symmetric makes CG does not converge
		const size_t DIM = size();
		for(size_t ci = 0; ci < DIM; ++ci) {
			Q_SetLen(&M_, ci+1, ptr[ci+1]-ptr[ci]);
			for(INT_TYPE nzi = ptr[ci]; nzi < ptr[ci+1]; ++nzi) {
//				if(idx[nzi] > ci) continue; // only upper triangular part
				Q_SetEntry(&M_, ci+1, nzi-ptr[ci], idx[nzi]+1, val[nzi]);
			}
		}
	}

	size_t size(void) const { return Q_GetDim(&M_); }

	::QMatrix &get(void) { return M_; }
	const ::QMatrix &get(void) const { return M_; }

	~QMatrix() {
		Q_Destr(&M_);
	}
private:
	mutable ::QMatrix M_;
};

class Vector
{
public:
	Vector(size_t DIM, const char *name) {
		V_Constr(&V_, const_cast<char *>(name), DIM, Normal, True);
	}
	template <typename T>
	void set(const T *v) {
		const size_t DIM = size();
		for(size_t i = 0; i < DIM; ++i)
			V_SetCmp(&V_, i+1, v[i]);
	}
	template <typename T>
	void get(T *v) const {
		for(size_t i = 0; i < size(); ++i)
			v[i] = V_GetCmp(&V_, i+1);
	}

	size_t size(void) const { return V_GetDim(&V_); }

	::Vector &get(void) { return V_; }
	const ::Vector &get(void) const { return V_; }


	~Vector() {
		V_Destr(&V_);
	}
private:
	mutable ::Vector V_;
};

};

size_t get_iter_type(const char *name[]);
IterProcType get_iter_proc(const char *name);

size_t get_precond_type(const char *name[]);
PrecondProcType get_precond_proc(const char *name);

#endif
