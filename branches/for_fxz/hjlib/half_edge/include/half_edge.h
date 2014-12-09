#ifndef HJ_HALF_EDGE_H_
#define HJ_HALF_EDGE_H_

//! @file half_edge.h
//! @brief the core file for half_edge data structure
//! @author Jin Huang
//! @date 201406-

//! @namespace hj::half_edge
//! @brief A concise and flexible half edge data structure

namespace hj { namespace half_edge {

//! @brief the minimal trival property
class default_interal_property
{
public:
  //! @brief property on vert
  struct vert {};
  //! @brief property on edge
  struct edge {};
  //! @brief property on face
  struct face {};
};

//! @brief wrap template argument for container
template <typename TMPL>
class mesh_tmpl
{
public:
  typedef TMPL mesh_t;
  mesh_t &operator()(void) {return *static_cast<mesh_t *>(this);}
  const mesh_t &operator()(void) const {return *static_cast<const mesh_t *>(this);}
};

//! @brief the half-edge data structure
template <template <typename> class ICon,
          typename IP=default_interal_property>
class half_edge_mesh_t : public mesh_tmpl<half_edge_mesh_t<ICon, IP> >
{
public: // types
  class vert_t;
  class edge_t;
  class face_t;

  typedef ICon<vert_t> verts_t;
  typedef ICon<edge_t> edges_t;
  typedef ICon<face_t> faces_t;

  //! @brief vert entity
  class vert_t : public IP::vert_t
  {
  public:
    typedef typename IP::vert_t prop_t;

    vert_t() {}
    vert_t(const prop_t &p):IP::vert_t(p) {}

    //@{
    //! @brief one of the edge pointing to this vertex
    const typename edges_t::const_iterator &edge() const { return edge_; }
    typename edges_t::const_iterator &edge() { return edge_; }
    //@}
  private:
    typename edges_t::const_iterator edge_;
  };

  //! @brief half_edge entity
  class edge_t : public IP::edge_t
  {
  public:
    typedef typename IP::edge_t prop_t;

    edge_t() {}
    edge_t(const prop_t &p):IP::edge_t(p) {}

    //@{
    //! @brief the vertex pointed by edge
    const typename verts_t::const_iterator &vert() const { return vert_; }
    typename verts_t::const_iterator &vert() { return vert_; }
    //@}

    //@{
    //! @brief the face neighboring to the half_edge
    const typename faces_t::const_iterator &face() const { return face_; }
    typename faces_t::const_iterator &face() { return face_; }
    //@}

    //@{
    //! @brief the next half_edge
    const typename edges_t::const_iterator &next() const { return next_; }
    typename edges_t::const_iterator &next() { return next_; }
    //@}

    //@{
    //! @brief the previous half_edge
    const typename edges_t::const_iterator &prev() const { return prev_; }
    typename edges_t::const_iterator &prev() { return prev_; }
    //@}

    //@{
    //! @brief the opposite half_edge
    const typename edges_t::const_iterator &oppo() const { return oppo_; }
    typename edges_t::const_iterator &oppo() { return oppo_; }
    //@}

  private:
    typename verts_t::const_iterator vert_;
    typename faces_t::const_iterator face_;
    typename edges_t::const_iterator next_, prev_, oppo_;
  };

  //! @brief face entity
  class face_t : public IP::face_t
  {
  public:
    typedef typename IP::face_t prop_t;

    face_t() {}
    face_t(const prop_t &p):IP::face_t(p) {}

    //@{
    //! @brief one of the edge in the face
    const typename edges_t::const_iterator &edge() const { return edge_; }
    typename edges_t::const_iterator &edge() { return edge_; }
    //@}
  private:
    typename edges_t::const_iterator edge_;
  };

public: // basic operations
  half_edge_mesh_t(){}

  //@{
  //! @brief assign from various possible types.
  template <typename M>
  half_edge_mesh_t(const M& m) {
    assign(*this, m);
  }
  template <typename M>
  const half_edge_mesh_t &operator = (const M& m) {
    assign(*this, m);
    return *this;
  }
  //@}

  //@{
  //! @brief get const container for query and traverse etc.
  const verts_t &verts() const { return verts_; }
  const edges_t &edges() const { return edges_; }
  const faces_t &faces() const { return faces_; }
  //@}

  //@{
  //! @brief add/del entity
  typename verts_t::const_iterator add(const vert_t &v) { return verts_.add(v); }
  typename edges_t::const_iterator add(const edge_t &e) { return edges_.add(e); }
  typename faces_t::const_iterator add(const face_t &f) { return faces_.add(f); }

  void del(const typename verts_t::const_iterator &v) { verts_.del(v); }
  void del(const typename edges_t::const_iterator &e) { edges_.del(e); }
  void del(const typename faces_t::const_iterator &f) { faces_.del(f); }
  //@}

  //@{
  //! @brief only way to change the entity pointed by the const_iterator
  vert_t &operator[](const typename verts_t::const_iterator &i) { return *const_cast<vert_t *>(i.operator->()); }
  edge_t &operator[](const typename edges_t::const_iterator &i) { return *const_cast<edge_t *>(i.operator->()); }
  face_t &operator[](const typename faces_t::const_iterator &i) { return *const_cast<face_t *>(i.operator->()); }
  //@}
private:
  verts_t verts_;
  edges_t edges_;
  faces_t faces_;
};

}}

#endif
