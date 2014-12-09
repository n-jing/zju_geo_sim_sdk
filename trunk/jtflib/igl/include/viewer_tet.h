#ifndef JTF_IGL_VIEWER_TET_H
#define JTF_IGL_VIEWER_TET_H

#include <igl/viewer/Viewer.h>
#include <Eigen/Core>

namespace igl
{
  class Viewer_tet : public igl::Viewer
  {
  public:
    int launch(std::string filename = "");

    bool load_tet_from_file(const char* mesh_file_name);
    bool save_tet_to_file(const char* mesh_file_name);

    Eigen::MatrixXi T; // Tetmesh #T x 4
  };
}
#endif // JTF_IGL_VIEWER_TET_H
