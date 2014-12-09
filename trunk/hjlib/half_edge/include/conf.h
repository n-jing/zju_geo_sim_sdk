#ifndef HJ_HALF_EDGE_CONF_H_
#define HJ_HALF_EDGE_CONF_H_

#ifdef __GNUG__
#define HJ_HALF_EDGE_DEPRECATED(msg) __attribute__((deprecated(msg)))
#endif

#ifdef _MSC_VER
#define HJ_HALF_EDGE_DEPRECATED(msg) __declspec(deprecated(msg))
#endif

#endif
