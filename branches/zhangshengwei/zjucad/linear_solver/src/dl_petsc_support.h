#ifndef _DL_PETSC_SUPPORT_H_
#define _DL_PETSC_SUPPORT_H_

#include "dl_petsc.h"

#if WIN32
#define LIB_PREFIX "lib"
#define LIB_SUFFIX ".so"
#else
#define LIB_PREFIX "lib"
#define LIB_SUFFIX ".so"
#endif
extern "C" {
  #include <gmodule.h>
}

#include <hjlib/sparse/sparse.h>

PETsc * ucreate_PETsc(const char * path = 0);
PETsc_CG * ucreate_PETsc_CG(const hj::sparse::csc<double, int32_t>& A, const char *pchar, const char * path=0);

#endif /*_DL_PETSC_SUPPORT_H_*/
