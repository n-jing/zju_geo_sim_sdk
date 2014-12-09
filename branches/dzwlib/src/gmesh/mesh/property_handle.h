#ifndef GMESH_MESH_PROPERTY_HANDLE_H_
#define GMESH_MESH_PROPERTY_HANDLE_H_

#include <vector>
#include "handle.h"

namespace gmesh {
// property handle
template <class T>
class prop_handle : public handle {

public:
  typedef T value;
  typedef T value_type;
  typedef std::vector<T> vector_type;
  typedef typename vector_type::reference reference;
  typedef typename vector_type::const_reference const_reference;

  explicit prop_handle(int index=-1) : handle(index) {}
};


template <class T>
class vprop_handle : public prop_handle<T> {
public:
  typedef T value;
  typedef T value_type;
  
  explicit vprop_handle(int index=-1) : prop_handle<T>(index) {}
  explicit vprop_handle(const prop_handle<T>& ph) : prop_handle<T> (ph) {}
};

template <class T>
class eprop_handle : public prop_handle<T> {
public:
  typedef T value;
  typedef T value_type;
  
  explicit eprop_handle(int index=-1) : prop_handle<T>(index) {}
  explicit eprop_handle(const prop_handle<T>& ph) : prop_handle<T> (ph) {}
};

template <class T>
class fprop_handle : public prop_handle<T> {
public:
  typedef T value;
  typedef T value_type;
  
  explicit fprop_handle(int index=-1) : prop_handle<T>(index) {}
  explicit fprop_handle(const prop_handle<T>& ph) : prop_handle<T> (ph) {}
};

template <class T>
class hprop_handle : public prop_handle<T> {
public:
  typedef T value;
  typedef T value_type;

  explicit hprop_handle(int index=-1) : prop_handle<T>(index) {}
  explicit hprop_handle(const prop_handle<T>& ph) : prop_handle<T> (ph) {}
};
}

#endif
