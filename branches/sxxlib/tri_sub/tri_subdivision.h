#ifndef TRI_SUBDIVISION_SXX_H
#define TRI_SUBDIVISION_SXX_H

#include <zjucad/matrix/matrix.h>

namespace sxx
{
  class Subdivision;
  int sub_edge(const zjucad::matrix::matrix<double> &points,
               const zjucad::matrix::matrix<size_t> &tri_faces,
               const std::vector<std::pair<size_t, size_t> > &edges,
               zjucad::matrix::matrix<double> &new_points,
               zjucad::matrix::matrix<size_t> &new_tri_faces,
               std::vector<size_t> &insert_points);

  int flip_edge(const zjucad::matrix::matrix<double> &points,
                const zjucad::matrix::matrix<size_t> &tri_faces,
                const std::vector<std::pair<size_t, size_t> > &edges,
                zjucad::matrix::matrix<double> &new_points,
                zjucad::matrix::matrix<size_t> &new_tri_faces,
                Subdivision *ptr = 0);
}

#endif
