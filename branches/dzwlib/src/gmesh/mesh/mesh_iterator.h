#ifndef GMESH_MESH_ITERATOR_H_
#define GMESH_MESH_ITERATOR_H_

namespace gmesh{

template<class mesh> class vv_iterator;
template<class mesh> class ve_iterator;
template<class mesh> class vf_iterator;  
template<class mesh> class vh_iterator;
template<class mesh> class fv_iterator;
template<class mesh> class fe_iterator;
template<class mesh> class ff_iterator;
template<class mesh> class fh_iterator;
template<class mesh> class const_vv_iterator;
template<class mesh> class const_ve_iterator;
template<class mesh> class const_vf_iterator;
template<class mesh> class const_vh_iterator;
template<class mesh> class const_fv_iterator;
template<class mesh> class const_fe_iterator;
template<class mesh> class const_ff_iterator;
template<class mesh> class const_fh_iterator;
/** \brief mesh item(vertex, edge, face) iterator 
 */
template<class mesh, class value_handle >
class mesh_iterator{  
 public:
  typedef value_handle value_type;
  typedef std::ptrdiff_t difference_type;
  typedef const value_type& reference;
  typedef const value_type* pointer;
  typedef std::bidirectional_iterator_tag iterator_category;
  typedef const mesh* mesh_ptr;
  typedef const mesh& mesh_ref;

 public:
  mesh_iterator(): mesh_(0) {}
  mesh_iterator(mesh_ref m, value_handle h) : mesh_(&m), handle_(h) {}

  reference operator*() const { return handle_; }
  pointer operator->() const { return &handle_; }
  value_handle get_handle() const { return handle_; }
  operator value_handle() const { return handle_; }

  bool operator==(const mesh_iterator& rhs) const {
    return mesh_ == rhs.mesh_ && handle_ == rhs.handle_;
  }
  bool operator!=(const mesh_iterator& rhs) const {
    return !operator==(rhs);
  }

  mesh_iterator& operator++() { handle_.increment(); }
  mesh_iterator& operator--() { handle_.decrement(); }
  
 private:
  mesh_ptr mesh_;
  value_handle handle_;
};
}

#include "vv_iterator.h"
#include "vf_iterator.h"
#include "ve_iterator.h"
#include "fv_iterator.h"
#include "fh_iterator.h"

#endif
