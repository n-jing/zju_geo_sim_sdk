#ifndef GMESH_MESH_FH_ITERATOR_H_
#define GMESH_MESH_FH_ITERATOR_H_

namespace gmesh {

template<class mesh> class fh_iterator;
template<class mesh> class const_fh_iterator;

template<class mesh>
class fh_iterator {
public:
  typedef typename mesh::halfedge_handle           halfedge_handle;
  typedef typename mesh::halfedge                  value_type;
  typedef typename mesh::halfedge_handle           value_handle;
  typedef typename mesh::halfedge&                 reference;
  typedef typename mesh::halfedge*                 pointer;
  typedef typename std::ptrdiff_t                  difference_type;
  typedef typename std::bidirectional_iterator_tag iterator_category;
  
  typedef mesh& mesh_ref;
  typedef mesh* mesh_ptr;

public:
  /// default constructor
  fh_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with mesh and face handle
  fh_iterator(mesh_ref m, typename mesh::face_handle fh) :
    mesh_(&m),
    start_heh_(m.get_halfedge_handle(fh)),
    curr_heh_(start_heh_),
    lap_counter_(0) {}
  /// construct with mesh and start halfedge
  fh_iterator(mesh_ref m, halfedge_handle heh) :
    mesh_(&m),
    start_heh_(heh),
    curr_heh_(heh),
    lap_counter_(0) {}
  /// copy constructor
  fh_iterator(const fh_iterator& rhs) :
    mesh_(rhs.mesh_),
    start_heh_(rhs.start_heh_),
    curr_heh_(rhs.curr_heh_),
    lap_counter_(0) {}
  /// assignment operator
  fh_iterator& operator=(const fh_iterator& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const fh_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_ &&
      curr_heh_ == rhs.curr_heh_ && lap_counter_ == rhs.lap_counter_;
  }
  bool operator!=(const fh_iterator& rhs) const {
    return !operator==(rhs);
  }

  fh_iterator& operator++() {
    assert(mesh_);
    curr_heh_ = mesh_->get_next_halfedge_handle(curr_heh_);
    if(curr_heh_ == start_heh_) lap_counter_++;
    return *this;
  }
  fh_iterator& operator--() {
    assert(mesh_);
    if(curr_heh_ == start_heh_) lap_counter_--;
    curr_heh_ = mesh_->get_prev_halfedge_handle(curr_heh_);
    return *this;
  }

  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  typename mesh::halfedge_handle handle() const {
    assert(mesh_);
    return curr_heh_;
  }

  operator typename mesh::halfedge_handle() const {
    assert(mesh_);
    return curr_heh_;
  }

  reference operator*() const {
    assert(mesh_);
    return mesh_->get_halfedge(curr_heh_);
  }

  pointer operator->() const {
    assert(mesh_);
    return &(mesh_->get_halfedge(curr_heh_));
  }

  operator bool() const {
    return curr_heh_.is_valid() && ((start_heh_ != curr_heh_) || (lap_counter_ == 0));
  }  

protected:
  mesh_ptr mesh_;
  halfedge_handle start_heh_, curr_heh_;
  int lap_counter_;
  
  friend class const_fh_iterator<mesh>;
};

template<class mesh>
class const_fh_iterator {
public:
  typedef typename mesh::halfedge_handle           halfedge_handle;
  typedef typename mesh::halfedge                  value_type;
  typedef typename mesh::halfedge_handle           value_handle;
  typedef const typename mesh::halfedge&           reference;
  typedef const typename mesh::halfedge*           pointer;
  typedef typename std::ptrdiff_t                  difference_type;
  typedef typename std::bidirectional_iterator_tag iterator_category;
  
  typedef const mesh& mesh_ref;
  typedef const mesh* mesh_ptr;

public:
  /// default constructor
  const_fh_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with mesh and face handle
  const_fh_iterator(mesh_ref m, typename mesh::face_handle fh) :
    mesh_(&m),
    start_heh_(m.get_halfedge_handle(fh)),
    curr_heh_(start_heh_),
    lap_counter_(0) {}
  /// construct with mesh and halfedge handle
  const_fh_iterator(mesh_ref m, halfedge_handle hh) :
    mesh_(&m),
    start_heh_(hh),
    curr_heh_(hh),
    lap_counter_(0) {}

  /// copy constructor
  const_fh_iterator(const const_fh_iterator& rhs):
  mesh_(rhs.mesh_),
  start_heh_(rhs.start_heh_),
  curr_heh_(rhs.curr_heh_),
  lap_counter_(rhs.lap_counter_) {}

  /// assignment operator
  const_fh_iterator& operator=(const const_fh_iterator& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  /// construct from non-const type
  const_fh_iterator(const fh_iterator<mesh>& rhs) :
    mesh_(rhs.mesh_),
    start_heh_(rhs.start_heh_),
    curr_heh_(rhs.curr_heh_),
    lap_counter_(rhs.lap_counter_) {}
  /// assign from non-const type
  const_fh_iterator& operator=(const fh_iterator<mesh>& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const const_fh_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_
      && curr_heh_ == rhs.curr_heh_ && lap_counter_ == rhs.lap_counter_;
  }
  bool operator!=(const const_fh_iterator& rhs) const {
    return !operator==(rhs);
  }

  /// next cw target
  const_fh_iterator& operator++() {
    assert(mesh_);
    curr_heh_ = mesh_->get_next_halfedge_handle(curr_heh_);
    if(curr_heh_ == start_heh_) lap_counter_++;
    return *this;
  }

  /// next ccw target
  const_fh_iterator& operator--() {
    assert(mesh_);
    if(curr_heh_ == start_heh_) lap_counter_--;
    curr_heh_ = mesh_->get_prev_halfedge_handle(curr_heh_);
    return *this;
  }
  
  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  /// return the handle of the current target
  typename mesh::halfedge_handle handle() const {
    assert(mesh_);
    return curr_heh_;
  }

  /// cast to the handle of the current target
  operator typename mesh::halfedge_handle() const {
    assert(mesh_);
    return curr_heh_;
  }

  /// return a reference to the current target
  reference operator*() const {
    assert(mesh_);
    return mesh_->get_halfedge(curr_heh_);
  }

  /// return a pointer to the current target
  pointer operator->() const {
    assert(mesh_);
    return *(mesh_->get_halfedge(curr_heh_));
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
