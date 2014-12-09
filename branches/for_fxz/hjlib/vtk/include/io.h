#ifndef HJ_VTK_IO_H_
#define HJ_VTK_IO_H_

#include <string>
#include <sstream>

#include <zjucad/matrix/matrix.h>

//! NOTE: ascii unstructured_grid only

namespace hj { namespace vtk {

inline std::string header(const char *desc) {
  return std::string("# vtk DataFile Version 2.0\n")
    + desc + "\nASCII\n";
}

inline std::string data_set(const char *type = "UNSTRUCTURED_GRID") {
  return std::string("DATASET ") + type + '\n';
}

template <typename OS, typename ITR>
OS& print(OS &os, const ITR beg, const ITR end, size_t num_per_row) {
  size_t c = 1;
  for(ITR i = beg; i != end; ++i, ++c) {
    os << *i;
    if(c == num_per_row) {
      os << '\n';
      c = 0;
    }
    else {
      os << ' ';
    }
  }
  return os;
}

using namespace zjucad::matrix;

template <typename E>
class vtk_points
{
public:
 vtk_points(const matrix_expression<E> &e)
    :e_(e){}
  const matrix_expression<E> &e_;
};

template <typename E>
vtk_points<E> points(const matrix_expression<E> &e) {
  return vtk_points<E>(e);
}

template <typename OS, typename E>
OS& operator << (OS &os, const vtk_points<E> &p) {
  os << "POINTS " << p.e_().size(2) << " FLOAT\n";
  return print(os, p.e_().begin(), p.e_().end(), p.e_().size(1));
}

template <typename E>
class vtk_cells
{
public:
 vtk_cells(const matrix_expression<E> &e)
    :e_(e){}
  const matrix_expression<E> &e_;
};

template <typename E>
vtk_cells<E> cells(const matrix_expression<E> &e) {
  return vtk_cells<E>(e);
}

template <typename OS, typename E>
OS& operator << (OS &os, const vtk_cells<E> &p) {
  const size_t points_per_cell = p.e_().size(1), num_cells = p.e_().size(2);
  os << "CELLS " << num_cells << ' ' << (points_per_cell+1)*num_cells << '\n';
  for(size_t i = 0; i < num_cells; ++i) {
    os << points_per_cell;
    for(size_t j = 0; j < points_per_cell; ++j)
      os << ' ' << p.e_()(j, i);
    os << '\n';
  }
  os << "CELL_TYPES " << num_cells << '\n';
  static const size_t guess_type[] = { // assume is simplex
    0, 1, 3, 5, 10
  };
  for(size_t i = 0; i < num_cells; ++i) {
    os << guess_type[points_per_cell] << '\n';
  }
  return os;
}

inline std::string point_data_begin(size_t np) {
  std::ostringstream oss;
  oss << "POINT_DATA " << np << '\n';
  return oss.str();
}

inline std::string cell_data_begin(size_t np) {
  std::ostringstream oss;
  oss << "CELL_DATA " << np << '\n';
  return oss.str();
}

template <typename E>
class vtk_data
{
public:
  vtk_data(const matrix_expression<E> &e, const char *name)
    :e_(e), name_(name) {}
  const matrix_expression<E> &e_;
  const std::string name_;
};

//! @param dim x num
template <typename E>
vtk_data<E> data(const matrix_expression<E>& e, const char *name) {
  return vtk_data<E>(e, name);
}

template <typename OS, typename E>
OS& operator << (OS &os, const vtk_data<E> &p) {
  if(p.e_().size(1) == 1) {
    os << "SCALARS " << p.name_ << " FLOAT 1\n"
       << "LOOKUP_TABLE default\n";
  }
  else {
    os << "VECTORS " << p.name_ << " FLOAT\n";
  }
  return print(os, p.e_().begin(), p.e_().end(), p.e_().size(1));
}

}}

#endif
