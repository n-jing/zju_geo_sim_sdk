#ifndef JTF_IGL_READTET_H
#define JTF_IGL_READTET_H
//#include <igl/igl_inline.h>

#ifndef IGL_NO_EIGEN
#include <Eigen/Core>
#endif
#endif

#include <string>
#include <vector>
#include <fstream>
#include <jtflib/mesh/io.h>
#include <igl/readOBJ.h>

namespace igl {
  template <typename Scalar, typename Index>
  bool readtet(
      const std::string tet_file_name,
      std::vector<std::vector<Scalar> > &V,
      std::vector<std::vector<Index> > &T)
  {
    using namespace std;
    using namespace zjucad::matrix;
    matrix<double> node;
    matrix<size_t> tet;
    if(jtf::mesh::tet_mesh_read_from_zjumat(tet_file_name.c_str(), &node, &tet))
      return false;
    V.resize(node.size(2));
    for(size_t ni = 0; ni < node.size(2); ++ni){
        V[ni].resize(node.size(1));
        for(size_t pi = 0; pi < node.size(1); ++pi){
            V[ni][pi] = node(pi,ni);
          }
      }
    T.resize(tet.size(2));
    for(size_t ti = 0; ti < tet.size(2); ++ti){
        T[ti].resize(tet.size(1));
        for(size_t pi = 0; pi < tet.size(1); ++pi){
            T[ti][pi] = tet(pi,ti);
          }
      }
    return true;
  }

#ifndef IGL_NO_EIGEN
  template <typename DerivedV, typename DerivedF>
  bool readtet(
      const std::string tet_file_name,
      Eigen::PlainObjectBase<DerivedV>& V,
      Eigen::PlainObjectBase<DerivedF>& T)
  {
    using namespace std;
    std::vector<std::vector<double> > node;
    std::vector<std::vector<int> > tet;
    bool success = igl::readtet(tet_file_name, node, tet);
    if(!success) return false;
    if(igl::list_to_matrix(node, V) == false) return false;
    if(igl::list_to_matrix(tet, T) == false)  return false;
    return true;
  }
#endif

}
