#ifndef GMESH_MESH_HALFEDGE_CONNECTIVITY_H_
#define GMESH_MESH_HALFEDGE_CONNECTIVITY_H_

#include <vector>
#include "connectivity.h"
#include "handle.h"

namespace gmesh {

class property_mananger;

class halfedge_connectivity : public connectivity{
public:
  typedef mesh_iterator<halfedge_connectivity, vert_handle> vert_iter;
  typedef mesh_iterator<halfedge_connectivity, edge_handle> edge_iter;
  typedef mesh_iterator<halfedge_connectivity, face_handle> face_iter;
  typedef mesh_iterator<halfedge_connectivity, halfedge_handle> halfedge_iter;

  typedef vert_iter const_vert_iter;
  typedef edge_iter const_edge_iter;
  typedef face_iter const_face_iter;
  typedef halfedge_iter const_halfedge_iter;

  typedef vv_iterator<halfedge_connectivity> vv_iter;
  typedef const_vv_iterator<halfedge_connectivity> const_vv_iter;
  typedef vf_iterator<gmesh::halfedge_connectivity> vf_iter;
  typedef const_vf_iterator<gmesh::halfedge_connectivity> const_vf_iter;
  typedef ve_iterator<halfedge_connectivity> ve_iter;
  typedef const_ve_iterator<halfedge_connectivity> const_ve_iter;
  typedef fv_iterator<halfedge_connectivity> fv_iter;
  typedef const_fv_iterator<halfedge_connectivity> const_fv_iter;
  typedef fh_iterator<halfedge_connectivity> fh_iter;
  typedef const_fh_iterator<halfedge_connectivity> const_fh_iter;

  //typedef gmesh::halfedge halfedge;
  //typedef gmesh::halfedge_handle halfedge_handle;

  vert_iter verts_begin() { return vert_iter(*this, vert_handle(0)); }
  vert_iter verts_end()   { return vert_iter(*this, vert_handle(pm_.n_verts())); }
  edge_iter edges_begin() { return edge_iter(*this, edge_handle(0)); }
  edge_iter edges_end()   { return edge_iter(*this, edge_handle(pm_.n_edges())); }
  face_iter faces_begin() { return face_iter(*this, face_handle(0)); }
  face_iter faces_end()   { return face_iter(*this, face_handle(pm_.n_faces())); }
  
  const_vert_iter verts_begin() const { return vert_iter(*this, vert_handle(0)); }
  const_vert_iter verts_end() const { return vert_iter(*this, vert_handle(pm_.n_verts())); }
  const_edge_iter edges_begin() const { return edge_iter(*this, edge_handle(0)); }
  const_edge_iter edges_end() const { return edge_iter(*this, edge_handle(pm_.n_edges())); }
  const_face_iter faces_begin() const { return face_iter(*this, face_handle(0)); }
  const_face_iter faces_end() const { return face_iter(*this, face_handle(pm_.n_faces())); }

  vv_iter vv_circulator(vert_handle vh) { return vv_iter(*this, vh); }
  const_vv_iter cvv_circulator(vert_handle vh) const { return const_vv_iter(*this, vh); }
  
private:
  typedef std::vector<halfedge> halfedge_container;  
  
public:
  halfedge_connectivity(property_mananger& pm) : connectivity(pm) {}
  virtual ~halfedge_connectivity() {}

  /** \brief add a new face and build the connectivity
   * Create a new face consisting of the vertices provided by vertex handle vector
   * These vertices have to be already added to the mesh
   * @param vhandles sorted list of vertex handles
   * (also defines order in which the vertices are added to the face)
   */
  virtual face_handle add_face(const std::vector<vert_handle>& vhandles);
  virtual edge_handle add_edge(vert_handle vh0, vert_handle vh1);

  virtual void delete_vert(vert_handle vh);
  virtual void delete_edge(edge_handle eh, bool is_del_iso_vert=true);
  virtual void delete_face(face_handle fh, bool is_del_iso_vert=true);
  virtual void garbage_collection();
  
  virtual bool is_boundary(vert_handle vh) const;
  virtual bool is_boundary(edge_handle eh) const;
  virtual bool is_boundary(face_handle fh) const;  
  
  virtual std::vector<face_handle> get_adj_face(vert_handle vh) const;
  virtual std::vector<face_handle> get_adj_face(edge_handle eh) const;
  virtual std::vector<face_handle> get_adj_face(face_handle fh) const;
  virtual std::vector<vert_handle> get_adj_vert(vert_handle vh) const;
  virtual std::vector<edge_handle> get_adj_edge(vert_handle vh) const;
  virtual std::vector<vert_handle> get_verts(edge_handle eh) const;
  virtual std::vector<vert_handle> get_verts(face_handle fh) const;
  virtual std::vector<edge_handle> get_edges(face_handle fh) const;
 
  bool is_boundary(halfedge_handle heh) const;
  halfedge_handle new_halfedge();
  halfedge_handle get_handle(const halfedge& he) const;
  halfedge& get_halfedge(halfedge_handle heh);
  const halfedge& get_halfedge(halfedge_handle heh) const;
  //void add_conn_base_property();
  bool is_valid_handle(halfedge_handle heh) const;
  bool is_valid_handle(vert_handle vh) const;
  bool is_valid_handle(edge_handle eh) const;
  bool is_valid_handle(face_handle fh) const;

  void set_boundary(halfedge_handle hh);

  // add by fxz
  bool is_flip_ok(edge_handle eh);
  void flip_edge(edge_handle eh);
  bool is_collapse_ok(halfedge_handle heh);
  void collapse_edge(halfedge_handle heh);
  
 public: 
  /** \brief find the halfedge hanlde
   * @param vh1 the start vertex handle
   * @param vh2 the end vertex handle
   * @return handle of halfedge
   */
  halfedge_handle find_halfedge(vert_handle vh1, vert_handle vh2) const;

  // -- get halfedge handle
  halfedge_handle get_prev_halfedge_handle(halfedge_handle heh) const;
  halfedge_handle get_next_halfedge_handle(halfedge_handle heh) const;
  halfedge_handle get_oppo_halfedge_handle(halfedge_handle heh) const;
  halfedge_handle get_halfedge_handle(vert_handle vh) const;
  halfedge_handle get_halfedge_handle(face_handle fh) const;
  halfedge_handle get_halfedge_handle(edge_handle eh, unsigned int i) const;

  halfedge_handle get_cw_rotated_halfedge_handle(halfedge_handle heh) const;
  halfedge_handle get_ccw_rotated_halfedge_handle(halfedge_handle heh) const;

  // -- get vertex handle
  vert_handle get_to_vert_handle(halfedge_handle heh) const;
  vert_handle get_from_vert_handle(halfedge_handle heh) const;

  vert_handle get_vert_handle(halfedge_handle heh) const; // vert
  face_handle get_face_handle(halfedge_handle heh) const; // face
  edge_handle get_edge_handle(halfedge_handle heh) const; // edge

  /// set vert/face/edge halfedge handle
  void set_halfedge_handle(vert_handle vh, halfedge_handle heh);
  void set_halfedge_handle(face_handle fh, halfedge_handle heh);
  void set_halfedge_handle(edge_handle eh, halfedge_handle heh0, halfedge_handle heh1);
  /// set halfedge's vert/face/prev_he/next_he/ handle
  void set_vert_handle(halfedge_handle heh, vert_handle vh);
  void set_edge_handle(halfedge_handle heh, edge_handle eh);
  void set_face_handle(halfedge_handle heh, face_handle fh);
  void set_prev_handle(halfedge_handle heh, halfedge_handle pre_heh);
  void set_next_handle(halfedge_handle heh, halfedge_handle nxt_heh);
  void set_oppo_handle(halfedge_handle heh, halfedge_handle opp_heh);
  
  /** \brief adjust a vertex' hh (halfedge handle) to be a boundary hh
   */ 
  void adjust_vert_halfedge_handle(vert_handle vh);

  status_info& get_status(vert_handle vh);
  const status_info& get_status(vert_handle vh) const;
  status_info& get_status(face_handle fh);
  const status_info& get_status(face_handle fh) const;
  status_info& get_status(edge_handle eh);
  const status_info& get_status(edge_handle eh) const;
  status_info& get_status(halfedge_handle hh);
  const status_info& get_status(halfedge_handle hh) const;

protected:
  void set_conn_base_property();
  const hprop_handle<status_info>& get_halfedge_status_property_handle() const { return hs_ph_; }
  hprop_handle<status_info>& get_halfedge_status_property_handle() { return hs_ph_; }

 private:
  hprop_handle<status_info> hs_ph_;
  
  struct edge_cache_info{
    halfedge_handle he_handle;
    bool is_new;
    bool needs_adjust;
  };
  std::vector<edge_cache_info> edge_cache_;
  // cache for set_next_halfedge and vertex's set_halfedge
  std::vector<std::pair<halfedge_handle, halfedge_handle> > next_cache_;
  unsigned int next_cache_count_;

  halfedge_container halfedges_;
};

inline void halfedge_connectivity::set_conn_base_property()
{
  pm_.add_property(get_halfedge_status_property_handle(), "h:status");  
}

inline halfedge_handle halfedge_connectivity::
get_handle(const halfedge& he) const {
  return halfedge_handle(&he - &halfedges_.front());
}

inline halfedge& halfedge_connectivity::
get_halfedge(halfedge_handle heh) {
  assert(is_valid_handle(heh));
  return halfedges_[heh.index()];
}

inline const halfedge& halfedge_connectivity::
get_halfedge(halfedge_handle heh) const {
  assert(is_valid_handle(heh));
  return halfedges_[heh.index()];
}


inline halfedge_handle halfedge_connectivity::
get_prev_halfedge_handle(halfedge_handle heh) const {
  return get_halfedge(heh).prev_he_handle_;
}

inline halfedge_handle halfedge_connectivity::
get_next_halfedge_handle(halfedge_handle heh) const {
  return get_halfedge(heh).next_he_handle_;
}

inline halfedge_handle halfedge_connectivity::
get_oppo_halfedge_handle(halfedge_handle heh) const {
  return halfedge_handle((heh.index() & 1) ?
                         heh.index()-1 : heh.index() + 1);
}

inline halfedge_handle halfedge_connectivity::
get_halfedge_handle(vert_handle vh) const {  
  return get_vert(vh).he_handle_;
}

inline halfedge_handle halfedge_connectivity::
get_halfedge_handle(face_handle fh) const {
  return get_face(fh).he_handle_;
}

inline halfedge_handle halfedge_connectivity::
get_halfedge_handle(edge_handle eh, unsigned int i) const {
  assert(i<=1); 
  return get_edge(eh).he_handles_[i];
}

inline vert_handle halfedge_connectivity::
get_to_vert_handle(halfedge_handle heh) const {
  return get_vert_handle(heh);
}

inline vert_handle halfedge_connectivity::
get_from_vert_handle(halfedge_handle heh) const {
  return get_to_vert_handle(get_oppo_halfedge_handle(heh));
}

inline vert_handle halfedge_connectivity::
get_vert_handle(halfedge_handle heh) const {
  return get_halfedge(heh).vert_handle_;
}

inline face_handle halfedge_connectivity::
get_face_handle(halfedge_handle heh) const {
  return get_halfedge(heh).face_handle_;
}

inline edge_handle halfedge_connectivity::
get_edge_handle(halfedge_handle heh) const {
  return get_halfedge(heh).edge_handle_;
}

inline halfedge_handle halfedge_connectivity::
get_cw_rotated_halfedge_handle(halfedge_handle heh) const{
  return get_next_halfedge_handle(get_oppo_halfedge_handle(heh));
}

inline halfedge_handle halfedge_connectivity::
get_ccw_rotated_halfedge_handle(halfedge_handle heh) const {
  return get_oppo_halfedge_handle(get_prev_halfedge_handle(heh));
}

inline void halfedge_connectivity::
set_halfedge_handle(vert_handle vh, halfedge_handle heh) {
  get_vert(vh).he_handle_ = heh;
}

inline void halfedge_connectivity::
set_halfedge_handle(face_handle fh, halfedge_handle heh) {
  get_face(fh).he_handle_ = heh;
}

inline void halfedge_connectivity::
set_halfedge_handle(edge_handle eh, halfedge_handle heh0, halfedge_handle heh1) {
  get_edge(eh).he_handles_[0] = heh0;
  get_edge(eh).he_handles_[1] = heh1;
}

inline void halfedge_connectivity::
set_vert_handle(halfedge_handle heh, vert_handle vh) {
  get_halfedge(heh).vert_handle_ = vh;
}

inline void halfedge_connectivity::
set_edge_handle(halfedge_handle heh, edge_handle eh) {
  get_halfedge(heh).edge_handle_ = eh;
}

inline void halfedge_connectivity::
set_face_handle(halfedge_handle heh, face_handle fh) {
  get_halfedge(heh).face_handle_ = fh;
}

inline void halfedge_connectivity::
set_prev_handle(halfedge_handle heh, halfedge_handle prev_heh) {
  get_halfedge(heh).prev_he_handle_ = prev_heh;
}

inline void halfedge_connectivity::
set_next_handle(halfedge_handle heh, halfedge_handle next_heh) {
  get_halfedge(heh).next_he_handle_ = next_heh;
}

inline bool halfedge_connectivity::is_boundary(halfedge_handle heh) const{
  return !get_face_handle(heh).is_valid();
}
inline bool halfedge_connectivity::is_boundary(edge_handle eh) const{
  return is_boundary(get_halfedge_handle(eh, 0)) || is_boundary(get_halfedge_handle(eh, 1));                     
}
inline bool halfedge_connectivity::is_boundary(face_handle fh) const{
  //FIXME:
  if (!fh.is_valid()) return true;
  std::vector<edge_handle> eh = get_edges(fh);
  for (int i=0; i<eh.size(); ++i)
    if (is_boundary(eh[i])) return true;
  
  return false;
}
inline bool halfedge_connectivity::is_boundary(vert_handle vh) const{
  //FIXME: is there a bug?
  halfedge_handle heh = get_halfedge_handle(vh);
  if (!heh.is_valid()) return true;
  do {
    if (is_boundary(heh)) return true;
    heh = get_cw_rotated_halfedge_handle(heh);
  } while (heh!=get_halfedge_handle(vh));

  return false;
}

inline void halfedge_connectivity::set_boundary(halfedge_handle hh) {
  return get_halfedge(hh).face_handle_.invalidate();
}
  
inline halfedge_handle halfedge_connectivity::new_halfedge() {
  halfedges_.push_back(halfedge());
  pm_.hprops_resize(halfedges_.size());
  return get_handle(halfedges_.back());
}

inline status_info& halfedge_connectivity::get_status(vert_handle vh) {
  return connectivity::get_status(vh);
}
inline const status_info& halfedge_connectivity::get_status(vert_handle vh) const {
  return connectivity::get_status(vh);
}
inline status_info& halfedge_connectivity::get_status(face_handle fh) {
  return connectivity::get_status(fh);
}
inline const status_info& halfedge_connectivity::get_status(face_handle fh) const {
  return connectivity::get_status(fh);
}
inline status_info& halfedge_connectivity::get_status(edge_handle eh) {
  return connectivity::get_status(eh);
}
inline const status_info& halfedge_connectivity::get_status(edge_handle eh) const {
  return connectivity::get_status(eh);
}

inline status_info& halfedge_connectivity::get_status(halfedge_handle hh) {
  return pm_.get_property(hs_ph_, hh);
}

inline const status_info& halfedge_connectivity::get_status(halfedge_handle hh) const {
  return pm_.get_property(hs_ph_, hh);
}

inline bool halfedge_connectivity::
is_valid_handle(halfedge_handle heh) const {
  return (heh.index() >= 0 && heh.index() < (int)halfedges_.size());
}

inline bool halfedge_connectivity::is_valid_handle(vert_handle vh) const
{
  return connectivity::is_valid_handle(vh);
}

inline bool halfedge_connectivity::is_valid_handle(edge_handle eh) const
{
  return connectivity::is_valid_handle(eh);
}

inline bool halfedge_connectivity::is_valid_handle(face_handle fh) const
{
  return connectivity::is_valid_handle(fh);
}

} //namespace

#endif
