#ifndef GMESH_MESH_PROPERTY_CONTAINER_H_
#define GMESH_MESH_PROPERTY_CONTAINER_H_

#include <vector>
#include <cassert>
#include <algorithm>
#include "property.h"
#include "property_handle.h"

namespace gmesh {

class property_container {
public:
  typedef std::vector<base_property*> properties;
  
  property_container();
  property_container(const property_container& pc);
  property_container& operator=(const property_container& rhs);
  ~property_container();  
  
  const properties& get_properties() const;
  properties& get_properties();
  const base_property& get_property(size_t index) const;
  base_property& get_property(size_t index);

  template<class T> const property<T>& get_property(prop_handle<T> h) const;
  template<class T> property<T>& get_property(prop_handle<T> h);

  template<class T> prop_handle<T> add_property(const T& prop, const std::string& name="<unknown>");
  template<class T> void remove_property(prop_handle<T> ph); // remove a property
  
  template<class T> prop_handle<T> get_handle(const T& prop, const std::string& name) const;

  void resize(size_t n) const;
  void reserve(size_t n) const;
  void swap(size_t i0, size_t i1) const;
  void clear();
  
  size_t size() const; // property size
  size_t n_elements() const; // element number of each property
  
 private:
  struct delete_functor {
    delete_functor(){}
    void operator() (base_property* p) const { if(p) delete p; p=NULL; }
  };
  struct resize_functor {
    resize_functor(size_t n) : n_(n) {}
    void operator() (base_property* p) const { if(p) p->resize(n_); }
    size_t n_;
  };
  struct reserve_functor {
    reserve_functor(size_t n) : n_(n) {}
    void operator() (base_property* p) const { if(p) p->reserve(n_); }
    size_t n_;
  };
  struct clear_functor {
    clear_functor() {}
    void operator() (base_property* p) const { if(p) p->clear(); }
  };
  struct swap_functor {
    swap_functor(size_t i0, size_t i1) : i0_(i0), i1_(i1) {}
    void operator() (base_property* p) const { if(p) p->swap(i0_, i1_);}
    size_t i0_, i1_;
  };
      
  properties properties_;
};

inline property_container::property_container() {}
inline property_container& property_container::operator=(const property_container& rhs) {
  // delete original properties
  std::for_each(properties_.begin(), properties_.end(), delete_functor());
  properties_ = rhs.properties_;
  properties::iterator p_it = properties_.begin(), p_end = properties_.end();
  for(; p_it != p_end; ++p_it)
    if(*p_it) *p_it = (*p_it)->clone();
  return *this;
}

inline property_container::property_container(const property_container& pc) {
  operator=(pc);
}

inline property_container::~property_container() {
  std::for_each(properties_.begin(), properties_.end(), delete_functor());
}

inline const property_container::properties& property_container::get_properties() const{
  return properties_;
}

inline property_container::properties& property_container::get_properties() { return properties_; }

inline const base_property& property_container::get_property(size_t index) const{
  assert(index < properties_.size());
  assert(properties_[index] != NULL);
  return *properties_[index];
}

inline base_property& property_container::get_property(size_t index) {
  assert(index < properties_.size());
  assert(properties_[index] != NULL);
  return *properties_[index];
}

template <class T>
inline const property<T>& property_container::get_property(prop_handle<T> h) const {
  assert(h.index() >= 0 && h.index() < (int) properties_.size());
  assert(properties_[h.index()] != NULL);
  property<T>* p = dynamic_cast<property<T>*> (properties_[h.index()]);
  assert(p != NULL);
  return *p;
}

template <class T>
inline property<T>& property_container::get_property(prop_handle<T> h){
  assert(h.index() >=0 && h.index() < (int) properties_.size());
  assert(properties_[h.index()] != NULL);
  property<T>* p = dynamic_cast<property<T>*> (properties_[h.index()]);
  assert(p != NULL);
  return *p;
}

template<class T>
inline prop_handle<T> property_container::add_property(const T&, const std::string& name){
  properties::iterator p_it = properties_.begin(), p_end = properties_.end();
  int index = 0;
  for(; p_it != p_end && *p_it != NULL; ++p_it, ++index);
  if(p_it == p_end) properties_.push_back(NULL);
  properties_[index] = new property<T>(name);
  return prop_handle<T>(index);
}

template<class T>
inline void property_container::remove_property(prop_handle<T> h) {
  assert(h.index() >=0 && h.index() < (int) properties_.size());
  assert(properties_[h.index()] != NULL);
  if(properties_[h.index()] != NULL) delete properties_[h.index()];
  properties_[h.index()] = NULL;
}

template<class T>
inline prop_handle<T> property_container::get_handle(const T&, const std::string& name) const {
  properties::const_iterator p_it = properties_.begin(), p_end = properties_.end();
  for(int index=0; p_it != p_end; ++p_it, ++index) {
    if((*p_it) != NULL && (*p_it)->name() == name &&
       dynamic_cast<property<T>*>(properties_[index] != NULL)) {
         return prop_handle<T>(index);
    }
  }
  
  return prop_handle<T>();
}

inline void property_container::resize(size_t n) const {
  std::for_each(properties_.begin(), properties_.end(), resize_functor(n));
}

// add by fxz
inline void property_container::clear() {
  std::for_each(properties_.begin(), properties_.end(), clear_functor());
}

inline size_t property_container::size() const {
  return properties_.size();
}

inline void property_container::reserve(size_t n) const {
  std::for_each(properties_.begin(), properties_.end(), reserve_functor(n));
}

 inline void property_container::swap(size_t i0, size_t i1) const {
   std::for_each(properties_.begin(), properties_.end(), swap_functor(i0, i1)); 
 }
 
inline size_t property_container::n_elements() const {
  if(properties_.size() == 0) return 0;
  return properties_[0]->n_elements();
}

}

#endif
