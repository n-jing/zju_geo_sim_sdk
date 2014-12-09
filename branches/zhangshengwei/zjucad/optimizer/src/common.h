#ifndef HJ_OPTIMIZER_COMMON_H_
#define HJ_OPTIMIZER_COMMON_H_

#include <cstdlib>
#include <stdint.h>

void dump(size_t iter_count, size_t dump_step, const char *dump_pref,
		  const double *x, int32_t nrow, int32_t ncol = 1);

#endif
