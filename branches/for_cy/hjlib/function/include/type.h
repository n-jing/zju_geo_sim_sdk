#ifndef HJ_FUNCTION_TYPE_H_
#define HJ_FUNCTION_TYPE_H_

#include <typeinfo>
#include <string>
#include <stdint.h>

#include "config.h"

namespace hj { namespace function {

template <typename T> inline char type2char(void) { return 0;}

template <> inline char type2char<float>(void) { return 's';}
template <> inline char type2char<double>(void) { return 'd';}
template <> inline char type2char<int32_t>(void) { return 'i';}
template <> inline char type2char<int64_t>(void) { return 'l';}

}}

#endif
