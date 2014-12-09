#ifndef GMESH_MESH_BASEPROPERTY_H_
#define GMESH_MESH_BASEPROPERTY_H_

#include <string>

namespace gmesh {

class base_property {

  static const size_t UnknownSize = size_t(-1);

 public:
  base_property(const std::string& name = "<unknown>") : name_(name){}
  base_property(const base_property& rhs) : name_(rhs.name_) {}
  virtual ~base_property() {}

 public: // synchronized array interface
  /// reserve memory for n elements
  virtual void reserve(size_t n) = 0;

  /// resize storage to hold n elements
  virtual void resize(size_t n) = 0;
  
  /// clear all elements and free memory
  virtual void clear() = 0;

  /// extend the number of elements by one
  virtual void push_back() = 0;

  virtual void swap(size_t i0, size_t i1) = 0;

  /// elements number 
  virtual size_t n_elements() const = 0;

  /// return a deep copy of self
  virtual base_property* clone() const = 0;

 public:
  const std::string& name() const { return name_; }

 private:
  std::string name_;
};

}

#endif
