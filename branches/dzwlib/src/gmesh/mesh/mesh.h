#ifndef GMESH_MESH_MESH_H_
#define GMESH_MESH_MESH_H_

#include "kernel.h"
#include "../common/types.h"

namespace gmesh {

template < class geometry_t, class connectivity_t>
class mesh : public kernel<geometry_t, connectivity_t>{
  
public:
  mesh(){}
  virtual ~mesh(){}

  vert_handle add_vertex(point p);
  face_handle add_face(const std::vector<vert_handle>& vhandles);
  face_handle add_face(vert_handle vh0, vert_handle vh1, vert_handle vh2);

  /// adjacent info api
 std::vector<face_handle> get_adj_face(vert_handle vh) const;
 std::vector<face_handle> get_adj_face(edge_handle eh) const;
 std::vector<face_handle> get_adj_face(face_handle fh) const;
 std::vector<vert_handle> get_adj_vert(vert_handle vh) const;
 std::vector<edge_handle> get_adj_edge(vert_handle eh) const;
 std::vector<vert_handle> get_verts(edge_handle eh) const;
 std::vector<vert_handle> get_verts(face_handle fh) const;
 //std::vector<vert_handle> get_verts(vert_handle vh) const;
 std::vector<edge_handle> get_edges(face_handle fh) const;

 point get_coord(vert_handle vh) const;
 void  set_coord(vert_handle vh, const point &p);

 ///
 size_t get_vert_num() const;
 size_t get_edge_num() const;
 size_t get_face_num() const;
  
  /// geometry api
  /// topology api
  
};

template<class geometry_t, class connectivity_t>
vert_handle mesh<geometry_t, connectivity_t>::add_vertex(point p){
  kernel<geometry_t, connectivity_t>::add_vert(p);
}

template<class geometry_t, class connectivity_t>
face_handle mesh<geometry_t, connectivity_t>::add_face(const std::vector<vert_handle>& vhandles) {
  kernel<geometry_t, connectivity_t>::add_face(vhandles);
}

template<class geometry_t, class connectivity_t>
face_handle mesh<geometry_t, connectivity_t>::add_face(vert_handle vh0, vert_handle vh1, vert_handle vh2) {
  std::vector<vert_handle> vhandles(3);
  vhandles[0] = vh0; vhandles[1] = vh1; vhandles[2] = vh2;
  kernel<geometry_t, connectivity_t>::add_face(vhandles);
}

template<class geometry_t, class connectivity_t>
std::vector<face_handle> mesh<geometry_t, connectivity_t>::get_adj_face(vert_handle vh) const {
  return kernel<geometry_t, connectivity_t>::get_adj_face(vh);
}

template<class geometry_t, class connectivity_t>
std::vector<face_handle> mesh<geometry_t, connectivity_t>::get_adj_face(edge_handle eh) const {
  return kernel<geometry_t, connectivity_t>::get_adj_face(eh);
}

template<class geometry_t, class connectivity_t>
std::vector<face_handle> mesh<geometry_t, connectivity_t>::get_adj_face(face_handle fh) const {
 return kernel<geometry_t, connectivity_t>::get_adj_face(fh);
}

template<class geometry_t, class connectivity_t>
std::vector<vert_handle> mesh<geometry_t, connectivity_t>::get_adj_vert(vert_handle vh) const {
 return kernel<geometry_t, connectivity_t>::get_adj_vert(vh);
}

template<class geometry_t, class connectivity_t>
std::vector<edge_handle> mesh<geometry_t, connectivity_t>::get_adj_edge(vert_handle vh) const {
 return kernel<geometry_t, connectivity_t>::get_adj_edge(vh);
}

/* template<class geometry_t, class connectivity_t> */
/* std::vector<vert_handle> mesh<geometry_t, connectivity_t>::get_verts(vert_handle vh) const { */
/*  return kernel<geometry_t, connectivity_t>::get_verts(vh); */
/* } */

template<class geometry_t, class connectivity_t>
std::vector<vert_handle> mesh<geometry_t, connectivity_t>::get_verts(edge_handle fh) const {
 return kernel<geometry_t, connectivity_t>::get_verts(fh);
}

template<class geometry_t, class connectivity_t>
std::vector<vert_handle> mesh<geometry_t, connectivity_t>::get_verts(face_handle fh) const {
 return kernel<geometry_t, connectivity_t>::get_verts(fh);
}

template<class geometry_t, class connectivity_t>
std::vector<edge_handle> mesh<geometry_t, connectivity_t>::get_edges(face_handle fh) const {
 return kernel<geometry_t, connectivity_t>::get_edges(fh);
}

template<class geometry_t, class connectivity_t>
point mesh<geometry_t, connectivity_t>::get_coord(vert_handle vh) const {
 return kernel<geometry_t, connectivity_t>::get_coord(vh);
}

template<class geometry_t, class connectivity_t>
    void mesh<geometry_t, connectivity_t>::set_coord(vert_handle vh, const point &p)
{
  return kernel<geometry_t, connectivity_t>::set_coord(vh, p);
}

template<class geometry_t, class connectivity_t>
size_t mesh<geometry_t, connectivity_t>::get_vert_num() const {
 return kernel<geometry_t, connectivity_t>::get_vert_num();
}

template<class geometry_t, class connectivity_t>
size_t mesh<geometry_t, connectivity_t>::get_edge_num() const {
 return kernel<geometry_t, connectivity_t>::get_edge_num();
}

template<class geometry_t, class connectivity_t>
size_t mesh<geometry_t, connectivity_t>::get_face_num() const {
 return kernel<geometry_t, connectivity_t>::get_face_num();
}

}

#endif
