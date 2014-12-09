#include "tri_subdivision.h"
#include "sxx_subdivision.h"

using namespace sxx;

namespace sxx
{
  int sub_edge(const zjucad::matrix::matrix<double> &points,
               const zjucad::matrix::matrix<size_t> &tri_faces,
               const std::vector<std::pair<size_t, size_t> > &edges,
               zjucad::matrix::matrix<double> &new_points,
               zjucad::matrix::matrix<size_t> &new_tri_faces,
               std::vector<size_t> &insert_points)
  {
    if(edges.size() == 0)
      {
        std::cerr << "the edge vector is empty " << std::endl;
        return 1;
      }
    Subdivision sub_mesh;
    sub_mesh.import_mesh(points, tri_faces);
    sub_mesh.subdived_edges(edges, insert_points);
    sub_mesh.export_mesh(new_points, new_tri_faces);
    return 0;
  }

  int flip_edge(const zjucad::matrix::matrix<double> &points,
                const zjucad::matrix::matrix<size_t> &tri_faces,
                const std::vector<std::pair<size_t, size_t> > &edges,
                zjucad::matrix::matrix<double> &new_points,
                zjucad::matrix::matrix<size_t> &new_tri_faces, Subdivision *ptr = 0)
  {
    if(edges.size() == 0)
      {
        std::cerr << "the edge vector is empty " << std::endl;
        return 1;
      }
    Subdivision sub_mesh;
    sub_mesh.import_mesh(points, tri_faces);
    sub_mesh.flip_edges(edges);
    sub_mesh.export_mesh(new_points, new_tri_faces);
    return 0;
  }
}
