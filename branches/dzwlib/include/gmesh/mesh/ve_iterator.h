#ifndef GMESH_MESH_VE_ITERATOR_H_
#define GMESH_MESH_VE_ITERATOR_H_

namespace gmesh {

template<class mesh> class ve_iterator;
template<class mesh> class const_ve_iterator;

template<class mesh>
class ve_iterator {
public:
  typedef typename mesh::halfedge_handle           halfedge_handle;
  typedef typename mesh::edge                      value_type;
  typedef typename mesh::edge_handle               value_handle;
  typedef typename mesh::edge&                     reference;
  typedef typename mesh::edge*                     pointer;
  typedef typename std::ptrdiff_t                  difference_type;
  typedef typename std::bidirectional_iterator_tag iterator_category;
  
  typedef mesh& mesh_ref;
  typedef mesh* mesh_ptr;

public:
  /// default constructor
  ve_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with mesh and vert handle
  ve_iterator(mesh_ref m, typename mesh::vert_handle vh) :
    mesh_(&m), start_heh_(m.get_halfedge_handle(vh)), curr_heh_(start_heh_), lap_counter_(0) {}
  /// construct with mesh and halfedge handle
  ve_iterator(mesh_ref m, halfedge_handle hh) :
    mesh_(&m), start_heh_(hh), curr_heh_(hh), lap_counter_(0) {}
  /// copy constructor
  ve_iterator(const ve_iterator& rhs) :
  mesh_(rhs.mesh_), start_heh_(rhs.start_heh_), curr_heh_(rhs.curr_heh_), lap_counter_(rhs.lap_counter_) {}
  /// assignment operator
  ve_iterator& operator=(const ve_iterator<mesh>& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const ve_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_ &&
      curr_heh_ == rhs.curr_heh_ && lap_counter_ == rhs.lap_counter_;
  }
  bool operator!=(const ve_iterator& rhs) const {
    return !operator==(rhs);
  } 

  ve_iterator& operator++() {
    assert(mesh_);
    curr_heh_ = mesh_->get_cw_rotated_halfedge_handle(curr_heh_);
    if(curr_heh_ == start_heh_) lap_counter_++;
    return *this;
  }

  ve_iterator& operator--() {
    assert(mesh_);
    if(curr_heh_ == start_heh_) lap_counter_--;
    curr_heh_ == mesh_->get_ccw_rotated_halfedge_handle(curr_heh_);
    return *this;    
  }

  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  typename mesh::edge_handle handle() const {
    assert(mesh_);
    return mesh_->get_edge_handle(curr_heh_);
  }

  operator typename mesh::edge_handle() const {
    assert(mesh_);
    return mesh_->get_edge_handle(curr_heh_);
  }

  reference operator*() const {
    assert(mesh_);
    return mesh_->get_edge(handle());
  }

  pointer operator->() const {
    assert(mesh_);
    return &(mesh_->get_edge(handle()));
  }

  operator bool() const {
    return curr_heh_.is_valid() && ((start_heh_ != curr_heh_) || (lap_counter_ == 0));
  }
  
protected:
  mesh_ptr mesh_;
  halfedge_handle start_heh_, curr_heh_;
  int lap_counter_;

  friend class const_ve_iterator<mesh>;
};
}

#endif
