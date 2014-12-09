#ifndef _DL_PETSC_SUPPORT_H_
#define _DL_PETSC_SUPPORT_H_

#include "dl_petsc.h"

PETsc * ucreate_PETsc(const char * path = 0);
PETsc_CG * ucreate_PETsc_CG(const double * val, const int32_t * idx,
			    const int32_t * ptr, const size_t nnz,
			    const size_t row,	const size_t col,
			    const char *pchar, const char * path=0);

#endif /*_DL_PETSC_SUPPORT_H_*/
