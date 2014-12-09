#ifndef GMESH_MESH_VV_ITERATOR_H_
#define GMESH_MESH_VV_ITERATOR_H_

namespace gmesh{

template<class mesh> class vv_iterator;
template<class mesh> class const_vv_iterator;

template<class mesh>
class vv_iterator{
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
  vv_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with the vertex handle
  vv_iterator(mesh_ref m, value_handle vh) :
      mesh_(&m), start_heh_(m.get_halfedge_handle(vh)),
      curr_heh_(start_heh_), lap_counter_(0) {}
  /// construct with the start halfedge handle
  vv_iterator(mesh_ref m, halfedge_handle hh) :
      mesh_(&m), start_heh_(hh),
      curr_heh_(hh), lap_counter_(0) {}

  vv_iterator(const vv_iterator& rhs) :
  mesh_(rhs.mesh_), start_heh_(rhs.start_heh_),
  curr_heh_(rhs.curr_heh_), lap_counter_(rhs.lap_counter_){}

  vv_iterator& operator=(const vv_iterator<mesh>& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const vv_iterator<mesh>& rhs) const {
    return (mesh_ == rhs.mesh_) && (start_heh_ == rhs.start_heh_) &&
        (curr_heh_ == rhs.curr_heh_) && (lap_counter_ == rhs.lap_counter_);
  }

  bool operator!=(const vv_iterator<mesh>& rhs) const {
    return !(operator==(rhs));
  }

  vv_iterator& operator++() {
    assert(mesh_);
    curr_heh_ = mesh_->get_ccw_halfedge_handle(curr_heh_);
    if(curr_heh_ == start_heh_) lap_counter_++;
    return *this;
  }

  vv_iterator& operator--() {
    assert(mesh_);
    if(curr_heh_ == start_heh_) lap_counter_--;
    curr_heh_ = mesh_->get_cw_halfedge_handle(curr_heh_);
    return *this;
  }

  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  vert_handle get_vert_handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }

  reference operator*() const {
    assert(mesh_);
    return mesh_->get_vert(get_vert_handle());
  }

  pointer operator->() const {
    assert(mesh_);
    return &(mesh_->get_vert(get_vert_handle()));
  }

  operator typename mesh::vert_handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }

  operator bool() const {
    return curr_heh_.is_valid() && ((start_heh_ != curr_heh_) || (lap_counter_ == 0));
  }
                                             
 private:
  mesh_ptr mesh_;
  halfedge_handle start_heh_; // the first halfedge handle
  halfedge_handle curr_heh_; // current halfedge handle
  int lap_counter_;

  friend class const_vv_iterator<mesh>;
};

template <class mesh>
class const_vv_iterator {
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

  const_vv_iterator(): mesh_(0), lap_counter_(0) {}
  const_vv_iterator(mesh_ref m, typename mesh::vert_handle vh):
      mesh_(&m), start_heh_(mesh_->get_halfedge_handle(vh)),
      curr_heh_(start_heh_), lap_counter_(0) {}

  const_vv_iterator(mesh_ref m, typename mesh::halfedge_handle hh) :
      mesh_(&m), start_heh_(hh), curr_heh_(hh), lap_counter_(0) {}

  const_vv_iterator(const const_vv_iterator& rhs):
  mesh_(rhs.mesh_), start_heh_(rhs.start_heh_),
  curr_heh_(rhs.curr_heh_), lap_counter_(rhs.lap_counter_) {}

  const_vv_iterator& operator=(const const_vv_iterator& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const const_vv_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_
        && curr_heh_ == rhs.curr_heh_ && lap_counter_ == rhs.lap_counter_;
  }

  bool operator!=(const const_vv_iterator& rhs) const {
    return !operator==(rhs);
  }

  reference operator*() const {
    assert(mesh_);
    return mesh_->get_vert(get_vert_handle());
  }

  pointer operator->() const {
    assert(mesh_);
    return &(mesh_->get_vert(get_vert_handle()));
  }

  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  vert_handle get_vert_handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }
  
  const_vv_iterator& operator++() {
    assert(mesh_);
    curr_heh_ = mesh_->get_cw_rotated_halfedge_handle(curr_heh_);
    if(curr_heh_ == start_heh_) lap_counter_++;
    return *this;
  }

  const_vv_iterator& operator--() {
    assert(mesh_);
    curr_heh_ = mesh_->get_ccw_rotated_halfedge_handle(curr_heh_);
    return *this;
  }
  
  operator typename mesh::vert_handle() const {
    assert(mesh_);
    return mesh_->get_to_vert_handle(curr_heh_);
  }

  operator bool() const {
    return curr_heh_.is_valid() && ((start_heh_ != curr_heh_) || (lap_counter_ == 0));
  }
  
 private:
  mesh_ptr mesh_;
  halfedge_handle start_heh_;
  halfedge_handle curr_heh_;
  int lap_counter_;
};

}

#endif
