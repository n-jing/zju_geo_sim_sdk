#ifndef GMESH_MESH_VF_ITERATOR_H_
#define GMESH_MESH_VF_ITERATOR_H_

#include <cassert>

namespace gmesh{

template<class mesh> class vf_iterator;
template<class mesh> class const_vf_iterator;

template<class mesh>  
class vf_iterator {
public:
  typedef typename mesh::halfedge_handle           halfedge_handle;
  typedef typename mesh::face                      value_type;
  typedef typename mesh::face_handle               value_handle;
  typedef typename mesh::face&                     reference;
  typedef typename mesh::face*                     pointer;
  typedef typename std::ptrdiff_t                  difference_type;
  typedef typename std::bidirectional_iterator_tag iterator_category;
  
  typedef mesh& mesh_ref;
  typedef mesh* mesh_ptr;

public:
  /// default constructor
  vf_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with mesh and vertex
  vf_iterator(mesh_ref m, typename mesh::vert_handle start_vh_) :
    mesh_(&m),
    start_heh_(m.get_halfedge_handle(start_vh_)),
    curr_heh_(start_heh_),
    lap_counter_(0) {}
  /// construct with mesh and start halfedge
  vf_iterator(mesh_ref m, halfedge_handle heh) :
    mesh_(&m),
    start_heh_(heh),
    curr_heh_(heh),
    lap_counter_(0) {}
  /// copy constructor
  vf_iterator(const vf_iterator& rhs) :
    mesh_(rhs.mesh_),
    start_heh_(rhs.start_heh_),
    curr_heh_(rhs.curr_heh_),
    lap_counter_(rhs.lap_counter_) {}
  /// asignment operator
  vf_iterator& operator=(const vf_iterator<mesh>& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  bool operator==(const vf_iterator& rhs) {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_ &&
      curr_heh_ == rhs.curr_heh_ && lap_counter_ == rhs.lap_counter_;
  }
  bool operator!=(const vf_iterator& rhs) {
    return !operator==(rhs);
  }

  /// pre_increment(next cw target)
  vf_iterator& operator++() {
    assert(mesh_);
    do {
      curr_heh_ = mesh_->get_cw_rotated_halfedge_handle(curr_heh_);
      if(curr_heh_ == start_heh_) lap_counter_ ++;
    }while ((*this) && (!handle().is_valid()));
    return *this;
  }
  /// pre-decrement(next ccw target)
  vf_iterator& operator--() {
    assert(mesh_);
    do {
      if(curr_heh_ == start_heh_) lap_counter_ --;
      curr_heh_ = mesh_->get_ccw_rotated_halfedge_handle(curr_heh_);
    } while ((*this) && (!handle().is_valid()));
  }
  
  halfedge_handle get_halfedge_handle() const { return curr_heh_; }

  /// return the handle of current target
  typename mesh::face_handle handle() const {
    assert(mesh_);
    return mesh_->get_face_handle(curr_heh_);
  }

  /// cast to the handle of the current target
  operator typename mesh::face_handle() const {
    assert(mesh_);
    return mesh_->get_face_handle(curr_heh_);
  }

  /// return a reference to the current target
  reference operator*() const {
    assert(mesh_);
    return mesh_->get_face(handle());
  }
  pointer operator->() const {
    assert(mesh_);
    return &(mesh_->get_face(handle()));
  }

  /** return whether the circulator is still valid.
      After one complete round around a vertex/face the ciruclator becomes
      invalid, i.e. this function will return \c false. Nevertheless you
      can continue circulating. This method just tells you whether you
      have completed the first round.
   */
  operator bool() const {
    return curr_heh_.is_valid() && ((start_heh_ != curr_heh_) || (lap_counter_ == 0));
  }
	       
  
protected:
  mesh_ptr mesh_;
  halfedge_handle start_heh_, curr_heh_;
  int lap_counter_;

  friend class const_vf_iterator<mesh>;
};

template<class mesh>
class const_vf_iterator {

public:
  typedef typename mesh::halfedge_handle           halfedge_handle;
  typedef typename mesh::face                      value_type;
  typedef typename mesh::face_handle               value_handle;
  typedef const typename mesh::face&                     reference;
  typedef const typename mesh::face*                     pointer;
  typedef typename std::ptrdiff_t                  difference_type;
  typedef typename std::bidirectional_iterator_tag iterator_category;
  
  typedef const mesh& mesh_ref;
  typedef const mesh* mesh_ptr;

public:
  /// default constructor
  const_vf_iterator() : mesh_(0), lap_counter_(0) {}
  /// construct with mesh and vertex handle
  const_vf_iterator(mesh_ref m, typename mesh::vert_handle start_vh) :
    mesh_(&m),
    start_heh_(m.get_halfedge_handle(start_vh)),
    curr_heh_(start_heh_),
    lap_counter_(0) {}
  /// copy constructor
  const_vf_iterator(const const_vf_iterator& rhs) :
  mesh_(rhs.mesh_),
  start_heh_(rhs.start_heh_),
  curr_heh_(rhs.curr_heh_),
  lap_counter_(rhs.lap_counter_) {}
  /// assignment operator
  const_vf_iterator& operator=(const const_vf_iterator<mesh>& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }

  /// construct from non-const vf_iterator
  const_vf_iterator(const vf_iterator<mesh>& rhs) :
    mesh_ (rhs.mesh_),
    start_heh_ (rhs.start_heh_),
    curr_heh_ (rhs.curr_heh_),
    lap_counter_ (rhs.lap_counter_) {}
  /// assign from non-const vf_iterator
  const_vf_iterator& operator=(const vf_iterator<mesh>& rhs) {
    mesh_ = rhs.mesh_;
    start_heh_ = rhs.start_heh_;
    curr_heh_ = rhs.curr_heh_;
    lap_counter_ = rhs.lap_counter_;
    return *this;
  }
  
  bool operator==(const const_vf_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && start_heh_ == rhs.start_heh_
      && curr_heh_ == rhs.curr_heh_  && lap_counter_ == rhs.lap_counter_;
  }
  bool operator!=(const const_vf_iterator& rhs) const {
    return !operator==(rhs);
  }

  /// pre-increment(next cw target)
  const_vf_iterator& operator++() {
    assert(mesh_);
    do {
      curr_heh_ = mesh_->get_cw_rotated_halfedge_handle(curr_heh_);
      if(curr_heh_ == start_heh_) lap_counter_ ++;      
    } while( (*this) && (!handle().is_valid()));
    return *this;
  }
  /// pre-decrement(next ccw target)
  const_vf_iterator& operator--() {
    assert(mesh_);
    do {
      if(curr_heh_ == start_heh_) lap_counter_--;
      curr_heh_ = mesh_->get_ccw_rotated_halfedge_handle(curr_heh_);
    } while( (*this) && (!handle().is_valid()));
    return *this;
  }

  halfedge_handle get_halfedge_handle() const {
    return curr_heh_;
  }

  typename mesh::face_handle handle() const {
    assert(mesh_);
    return mesh_->get_face_handle(curr_heh_);
  }

  operator typename mesh::face_handle() const {
    assert(mesh_);
    return mesh_->get_face_handle(curr_heh_);
  }

  reference operator*() const {
    assert(mesh_);
    return mesh_->get_face(handle());
  }
  pointer operator->() const {
    assert(mesh_);
    return mesh_->get_face(handle());
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
