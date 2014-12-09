#ifndef HJ_HALF_EDGE_IO_H_
#define HJ_HALF_EDGE_IO_H_

#include "builder.h"
#include <jtflib/mesh/io.h>

template <typename TMPL>
bool read_mesh(TMPL & m, const char* _filename);


template <typename TMPL>
bool write_mesh(TMPL & m, const char* _filename);

#include "IO.hh"

#endif
