#ifndef GMESH_MESH_FV_ITERATOR_H_
#define GMESH_MESH_FV_ITERATOR_H_

namespace gmesh {

template<class mesh> class fv_iterator;
template<class mesh> class const_fv_iterator;

template<class mesh>
class fv_iterator {
public:
  typedef typename mesh::halfedge_handle           halfedge_handle;
  typedef typename mesh::vert                      value_type;
  typedef typename mesh::vert_handle               value_handle;
  typedef typename mesh::vert&                     reference;
  typedef typename mesh::vert*                     pointer;
  typedef typename std::ptrdiff_t                  difference_type;
  typedef typename std::bidirectional_iterator_tag iterator_category;
  
  typedef mesh& mesh_ref;
  typedef mesh* mesh_ptr;

public:
  /// default constructor
  fv_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with mesh and face handle
  fv_iterator(mesh_ref m, typename mesh::face_handle fh) :
    mesh_(&m),
    start_heh_(m.get_halfedge_handle(fh)),
    curr_heh_(start_heh_),
    lap_counter_(0) {}
  /// construct with mesh and start halfedge
  fv_iterator(mesh_ref m, halfedge_handle heh) :
    mesh_(&m),
    start_heh_(heh),
    curr_heh_(heh),
    lap_counter_(0) {}
  /// copy constructor
  fv_iterator(const fv_iterator& rhs) :
  mesh_(rhs.mesh_),
  start_heh_(rhs.start_heh_),
  curr_heh_(rhs.curr_heh_),
  lap_counter_(0) {}
  /// assignment operator
  fv_iterator& operator=(const fv_iterator& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const fv_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_ &&
      curr_heh_ == rhs.curr_heh_ && lap_counter_ == rhs.lap_counter_;
  }
  bool operator!=(const fv_iterator& rhs) const {
    return !operator==(rhs);
  }

  fv_iterator& operator++() {
    assert(mesh_);
    curr_heh_ = mesh_->get_next_halfedge_handle(curr_heh_);
    if(curr_heh_ == start_heh_) lap_counter_++;
    return *this;
  }
  fv_iterator& operator--() {
    assert(mesh_);
    if(curr_heh_ == start_heh_) lap_counter_--;
    curr_heh_ = mesh_->get_prev_halfedge_handle(curr_heh_);
    return *this;
  }

  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  typename mesh::vert_handle handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }

  operator typename mesh::vert_handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }

  reference operator*() const {
    assert(mesh_);
    return mesh_->get_vert(handle());
  }

  pointer operator->() const {
    assert(mesh_);
    return &(mesh_->get_vert(handle()));
  }

  operator bool() const {
    return curr_heh_.is_valid() && ((start_heh_ != curr_heh_) || (lap_counter_ == 0));
  }  

protected:
  mesh_ptr mesh_;
  halfedge_handle start_heh_, curr_heh_;
  int lap_counter_;
  
  friend class const_fv_iterator<mesh>;
};

template<class mesh>
class const_fv_iterator {
public:
  typedef typename mesh::halfedge_handle           halfedge_handle;
  typedef typename mesh::vert                      value_type;
  typedef typename mesh::vert_handle               value_handle;
  typedef const typename mesh::vert&               reference;
  typedef const typename mesh::vert*               pointer;
  typedef typename std::ptrdiff_t                  difference_type;
  typedef typename std::bidirectional_iterator_tag iterator_category;
  
  typedef const mesh& mesh_ref;
  typedef const mesh* mesh_ptr;

public:
  /// default constructor
  const_fv_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with mesh and face handle
  const_fv_iterator(mesh_ref m, typename mesh::face_handle fh) :
    mesh_(&m),
    start_heh_(m.get_halfedge_handle(fh)),
    curr_heh_(start_heh_),
    lap_counter_(0) {}
  /// construct with mesh and halfedge handle
  const_fv_iterator(mesh_ref m, halfedge_handle hh) :
    mesh_(&m),
    start_heh_(hh),
    curr_heh_(hh),
    lap_counter_(0) {}

  /// copy constructor
  const_fv_iterator(const const_fv_iterator& rhs):
  mesh_(rhs.mesh_),
  start_heh_(rhs.start_heh_),
  curr_heh_(rhs.curr_heh_),
  lap_counter_(rhs.lap_counter_) {}

  /// assignment operator
  const_fv_iterator& operator=(const const_fv_iterator& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  /// construct from non-const type
  const_fv_iterator(const fv_iterator<mesh>& rhs) :
    mesh_(rhs.mesh_),
    start_heh_(rhs.start_heh_),
    curr_heh_(rhs.curr_heh_),
    lap_counter_(rhs.lap_counter_) {}
  /// assign from non-const type
  const_fv_iterator& operator=(const fv_iterator<mesh>& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const const_fv_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_
      && curr_heh_ == rhs.curr_heh_ && lap_counter_ == rhs.lap_counter_;
  }
  bool operator!=(const const_fv_iterator& rhs) const {
    return !operator==(rhs);
  }

  /// next cw target
  const_fv_iterator& operator++() {
    assert(mesh_);
    curr_heh_ = mesh_->get_next_halfedge_handle(curr_heh_);
    if(curr_heh_ == start_heh_) lap_counter_++;
    return *this;
  }

  /// next ccw target
  const_fv_iterator& operator--() {
    assert(mesh_);
    if(curr_heh_ == start_heh_) lap_counter_--;
    curr_heh_ = mesh_->get_prev_halfedge_handle(curr_heh_);
    return *this;
  }
  
  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  /// return the handle of the current target
  typename mesh::vert_handle handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }

  /// cast to the handle of the current target
  operator typename mesh::vert_handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }

  /// return a reference to the current target
  reference operator*() const {
    assert(mesh_);
    return mesh_->get_vert(handle());
  }

  /// return a pointer to the current target
  pointer operator->() const {
    assert(mesh_);
    return *(mesh_->get_vert(handle()));
  }

  operator bool() const {
    return curr_heh_.is_valid() && ((start_heh_ != curr_heh_) || (lap_counter_ == 0));
  }

protected:
  mesh_ptr mesh_;
  halfedge_handle start_heh_, curr_heh_;
  int lap_counter_;
};
  
}

#endif
