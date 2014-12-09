#ifndef GMESH_MESH_MESH_ITEM_H_
#define GMESH_MESH_MESH_ITEM_H_

#include "handle.h"

namespace gmesh {
struct vert {
  halfedge_handle he_handle_;
};

struct face {
  halfedge_handle he_handle_;
};

struct edge {
  halfedge_handle he_handles_[2];
};

struct halfedge {
  face_handle face_handle_;
  edge_handle edge_handle_;
  vert_handle vert_handle_;
  halfedge_handle prev_he_handle_;
  halfedge_handle next_he_handle_;
  halfedge_handle oppo_he_handle_;
};
}

#endif
