#ifndef GMESH_MESH_PROPERTY_H_
#define GMESH_MESH_PROPERTY_H_

#include "base_property.h"

namespace gmesh {

template <class T>
class property : public base_property {  

public:
  typedef T                                     value_type;
  typedef std::vector<T>                        vector_type;
  typedef typename vector_type::reference       reference;
  typedef typename vector_type::const_reference const_reference;
  
public:
  property(const std::string& name = "<unknown>") : base_property(name) {}
  property(const property<T>& rhs) : base_property( rhs), data_(rhs.data_) {}

public:
  virtual void reserve(size_t n) { data_.reserve(n); }
  virtual void resize(size_t n) { data_.resize(n); }
  virtual void clear() { data_.clear(); vector_type().swap(data_); }
  virtual void push_back() { data_.push_back(T()); }
  virtual void swap(size_t i0, size_t i1) { std::swap(data_[i0], data_[i1]); }

public:
  virtual size_t n_elements() const { return data_.size(); }


public:
  const T* data() const {
    if( data_.empty() )
      return 0;

    return &data_[0];
  }

  reference operator[] (int idx ) {
    assert( size_t(idx) < data_.size() );
    return data_[idx];
  }

  const_reference operator[] (int idx) const {
    assert( size_t(idx) < data_.size() );
    return data_[idx];
  }

  property<T>* clone() const {
    property<T>* p = new property<T> (*this);
    return p;
  }

private:
  vector_type data_;

};

}

#endif
