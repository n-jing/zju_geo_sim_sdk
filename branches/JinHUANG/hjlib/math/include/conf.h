#ifndef HJ_MATH_CONF_H_
#define HJ_MATH_CONF_H_

#ifdef WIN32
#  ifdef hj_math_EXPORTS
#    define HJ_MATH_API __declspec(dllexport)
#  else
#    define HJ_MATH_API
#  endif
#else
#  define HJ_MATH_API
#endif

#endif
