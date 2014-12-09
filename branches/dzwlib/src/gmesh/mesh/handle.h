#ifndef GMESH_MESH_HANDLE_H_
#define GMESH_MESH_HANDLE_H_

namespace gmesh {
/*! Basic handle class, inherted by other handle classes
*/
class handle {
  
 public:
  explicit handle(int index=-1) : index_(index) {}

  bool operator ==(const handle& rhs) const { return index_ == rhs.index_; }
  bool operator !=(const handle& rhs) const { return index_ != rhs.index_; }
  int index() const { return index_; }
  bool is_valid() const { return index_ != -1; }
  void reset() { index_ = -1; }
  void invalidate() { index_ = -1; }

  // this is to be used only by iterator
  void increment() { ++index_; }
  void decrement() { --index_; }

 private:
  int index_;
};  

class vert_handle : public handle {
 public:
  explicit vert_handle(int index=-1) : handle(index) {}
};

class face_handle : public handle {
 public:
  explicit face_handle(int index=-1) : handle(index) {}
};

class edge_handle : public handle {
 public:
  explicit edge_handle(int index=-1) : handle(index) {}
};

class halfedge_handle : public handle {
public:
  explicit halfedge_handle(int index=-1) : handle(index) {}
};
}

#endif
