#ifndef GMESH_MESH_CONNECTIVITY_H_
#define GMESH_MESH_CONNECTIVITY_H_

#include "handle.h"
#include "mesh_item.h"
#include "mesh_iterator.h"
#include "property_mananger.h"
#include "property_handle.h"
#include "item_status.h"

#include <vector>
#include <cassert>

namespace gmesh {

class connectivity{
public:
  typedef gmesh::vert_handle vert_handle;
  typedef gmesh::edge_handle edge_handle;
  typedef gmesh::face_handle face_handle;
  typedef gmesh::halfedge_handle halfedge_handle;
  
  typedef gmesh::vert vert;
  typedef gmesh::edge edge;
  typedef gmesh::face face;
  typedef gmesh::halfedge halfedge;
  
private:
  typedef std::vector<vert> vert_container;
  typedef std::vector<edge> edge_container;
  typedef std::vector<face> face_container;
  
public:
  connectivity(property_mananger& pm) : pm_(pm) {}
  virtual ~connectivity() {}
  
  virtual face_handle add_face(const std::vector<vert_handle>& vhandles) = 0;
  virtual edge_handle add_edge(vert_handle vh0, vert_handle vh1) = 0;

  virtual void delete_vert(vert_handle vh) = 0;
  virtual void delete_edge(edge_handle eh, bool is_del_iso_vert) = 0;
  virtual void delete_face(face_handle fh, bool is_del_iso_vert) = 0;
  virtual void garbage_collection() = 0; /// collect deleted items

  virtual vert_handle new_vert();
  virtual edge_handle new_edge();
  virtual face_handle new_face();
  

  virtual bool is_boundary(vert_handle vh) const = 0;
  virtual bool is_boundary(edge_handle eh) const = 0;
  virtual bool is_boundary(face_handle fh) const = 0;
  
  virtual vert_handle get_handle(const vert& v) const;
  virtual edge_handle get_handle(const edge& e) const;
  virtual face_handle get_handle(const face& f) const;

  const vert& get_vert(vert_handle vh) const;
  vert& get_vert(vert_handle vh);
  const edge& get_edge(edge_handle eh) const;
  edge& get_edge(edge_handle eh);
  const face& get_face(face_handle fh) const;
  face& get_face(face_handle fh);  
  
  virtual vert_handle get_vert_handle(unsigned int i) const;
  virtual edge_handle get_edge_handle(unsigned int i) const;
  virtual face_handle get_face_handle(unsigned int i) const;

  virtual std::vector<face_handle> get_adj_face(vert_handle vh) const = 0;
  virtual std::vector<face_handle> get_adj_face(edge_handle eh) const = 0;
  virtual std::vector<face_handle> get_adj_face(face_handle fh) const = 0;
  virtual std::vector<vert_handle> get_adj_vert(vert_handle vh) const = 0;
  virtual std::vector<edge_handle> get_adj_edge(vert_handle vh) const = 0;
  virtual std::vector<vert_handle> get_verts(face_handle fh) const = 0;
  virtual std::vector<vert_handle> get_verts(edge_handle eh) const = 0;
  virtual std::vector<edge_handle> get_edges(face_handle fh) const = 0;

  virtual bool is_valid_handle(vert_handle vh) const;
  virtual bool is_valid_handle(edge_handle eh) const;
  virtual bool is_valid_handle(face_handle fh) const;

  size_t get_vert_num() const;
  size_t get_edge_num() const;
  size_t get_face_num() const;

  /// get status info
  virtual status_info& get_status(vert_handle vh);
  virtual const status_info& get_status(vert_handle vh) const;
  virtual status_info& get_status(face_handle fh);
  virtual const status_info& get_status(face_handle fh) const;
  virtual status_info& get_status(edge_handle eh);
  virtual const status_info& get_status(edge_handle eh) const;
  //virtual status_info& get_status(halfedge_handle hh);
  //virtual const status_info& get_status(halfedge_handle hh) const ;

  /// get status properties
  virtual vprop_handle<status_info>& get_vert_status_property_handle() { return vs_ph_; }
  virtual const vprop_handle<status_info>& get_vert_status_property_handle() const { return vs_ph_; }
  virtual fprop_handle<status_info>& get_face_status_property_handle() { return fs_ph_; }
  virtual const fprop_handle<status_info>& get_face_status_property_handle() const { return fs_ph_; }
  virtual eprop_handle<status_info>& get_edge_status_property_handle() { return es_ph_; }
  virtual const eprop_handle<status_info>& get_edge_status_property_handle() const { return es_ph_; }

protected:
  vert_container verts_;
  edge_container edges_;
  face_container faces_;

  /// status property
  vprop_handle<status_info> vs_ph_;
  fprop_handle<status_info> fs_ph_;
  eprop_handle<status_info> es_ph_;

protected:
  property_mananger& pm_;
};

inline vert_handle connectivity::get_handle(const vert& v) const {
  return vert_handle(&v - &verts_.front());
}

inline edge_handle connectivity::get_handle(const edge& e) const {
  return edge_handle(&e - &edges_.front());
}

inline face_handle connectivity::get_handle(const face& f) const {
  return face_handle(&f - &faces_.front());
}

inline const vert& connectivity::get_vert(vert_handle vh) const {
  assert(is_valid_handle(vh));
  return verts_[vh.index()];
}

inline vert& connectivity::get_vert(vert_handle vh) {
  assert(is_valid_handle(vh));
  return verts_[vh.index()];
}

inline const edge& connectivity::get_edge(edge_handle eh) const {
  assert(is_valid_handle(eh));
  return edges_[eh.index()];
}

inline edge& connectivity::get_edge(edge_handle eh) {
  assert(is_valid_handle(eh));
  return edges_[eh.index()];
}

inline const face& connectivity::get_face(face_handle fh) const {
  assert(is_valid_handle(fh));
  return faces_[fh.index()];
}

inline face& connectivity::get_face(face_handle fh) {
  assert(is_valid_handle(fh));
  return faces_[fh.index()];
}

inline vert_handle connectivity::get_vert_handle(unsigned int i) const {
  assert(i < verts_.size());
  return get_handle(verts_[i]);
}

inline edge_handle connectivity::get_edge_handle(unsigned int i) const {
  assert(i < edges_.size());
  return get_handle(edges_[i]);
}

inline face_handle connectivity::get_face_handle(unsigned int i) const {
  assert(i < faces_.size());
  return get_handle(faces_[i]);
}

inline bool connectivity::is_valid_handle(vert_handle vh) const {
  return vh.index() >= 0 && vh.index() < verts_.size();
}

inline bool connectivity::is_valid_handle(edge_handle eh) const {
  return eh.index() >= 0 && eh.index() < edges_.size();
}

inline bool connectivity::is_valid_handle(face_handle fh) const {
  return fh.index() >= 0 && fh.index() < faces_.size();
}

inline size_t connectivity::get_vert_num() const {
  return verts_.size();
}

inline size_t connectivity::get_edge_num() const {
  return edges_.size();
}

inline size_t connectivity::get_face_num() const {
  return faces_.size();
}

inline vert_handle connectivity::new_vert() {
  verts_.push_back(vert());
  pm_.vprops_resize(verts_.size());
  return get_handle(verts_.back());
}

inline edge_handle connectivity::new_edge() {
  edges_.push_back(edge());
  pm_.eprops_resize(edges_.size());
  return get_handle(edges_.back());
}

inline face_handle connectivity::new_face() {
  faces_.push_back(face());
  pm_.fprops_resize(faces_.size());
  return get_handle(faces_.back());
}

inline status_info& connectivity::get_status(vert_handle vh) {
  return pm_.get_property(vs_ph_, vh);
}
inline const status_info& connectivity::get_status(vert_handle vh) const {
  return pm_.get_property(vs_ph_, vh);
}
inline status_info& connectivity::get_status(face_handle fh) {
  return pm_.get_property(fs_ph_, fh);
}
inline const status_info& connectivity::get_status(face_handle fh) const {
  return pm_.get_property(fs_ph_, fh);
}
inline status_info& connectivity::get_status(edge_handle eh) {
  return pm_.get_property(es_ph_, eh);
}
inline const status_info& connectivity::get_status(edge_handle eh) const {
  return pm_.get_property(es_ph_, eh);
}
}

#endif 
