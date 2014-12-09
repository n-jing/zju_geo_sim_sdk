#ifndef GMESH_MESH_PROPERTY_MANANGER_H_
#define GMESH_MESH_PROPERTY_MANANGER_H_

#include "property_container.h"

#include <string>
#include <vector>

namespace gmesh {
class property_mananger{
public:
  typedef unsigned int uint;
  
public:
  property_mananger() {}
  virtual ~property_mananger(){}


  /// \name add a property to a mesh item
  /** \brief add a vertex property
   * \param ph A propert handle defining the data type to bind to mesh.
   *           On success the handle is valid else unvalid
   * \param name Optional name of property
   */
  template<class T>
  void add_property(vprop_handle<T>& ph, const std::string name="<vprop>"){
    ph = vprop_handle<T>(vprops_.add_property(T(), name));
    vprops_.resize(vprops_.n_elements());
  }
  /// \name remove a property from a mesh item
  template<class T>
  void remove_property(vprop_handle<T>& ph) {
    if ( ph.is_valid()) {
      vprops_.remove_property(ph);
      ph.reset();
    }
  }

  template<class T>
  void add_property(fprop_handle<T>& ph, const std::string name="<fprop>"){
    ph = fprop_handle<T>(fprops_.add_property(T(), name));
    fprops_.resize(fprops_.n_elements());
  }
  template<class T>
  void remove_property(fprop_handle<T>& ph) {
    if ( ph.is_valid()) {
      fprops_.remove_property(ph);
      ph.reset();
    }
  }

  template<class T>
  void add_property(eprop_handle<T>& ph, const std::string name="<fprop>"){
    ph = eprop_handle<T>(eprops_.add_property(T(), name));
    eprops_.resize(eprops_.n_elements());
  }
  template<class T>
  void remove_property(eprop_handle<T>& ph) {
    if (ph.is_valid()) {
      eprops_.remove_property(ph);
      ph.reset();
    }
  }

  template<class T>
  void add_property(hprop_handle<T>& ph, const std::string name="<hprop>"){
    ph = hprop_handle<T>(hprops_.add_property(T(), name));
    hprops_.resize(hprops_.n_elements());
  }
  template<class T>
  void remove_property(hprop_handle<T>& ph) {
    if (ph.is_valid()) {
      hprops_.remove_property(ph);
      ph.reset();
    }
  }
  

public: // get handle from name
  template<class T>
  bool get_property_handle(vprop_handle<T>& ph, const std::string& name) const {
    return (ph = vprop_handle<T>(vprops_.get_handle(T(), name))).is_valid();
  }

  template<class T>
  bool get_property_handle(fprop_handle<T>& fh, const std::string& name) const {
    return (fh = fprop_handle<T>(fprops_.get_handle(T(), name))).is_valid();
  }

  template<class T>
  bool get_property_handle(eprop_handle<T>& eh, const std::string& name) const {
    return (eh = eprop_handle<T>(eprops_.get_handle(T(), name))).is_valid();
  }

  template<class T>
  bool get_property_handle(hprop_handle<T>& hh, const std::string& name) const {
    return (hh = hprop_handle<T>(hprops_.get_handle(T(), name))).is_valid();
  }
  
public: // access mesh property
  template <class T>
  property<T>& get_property(vprop_handle<T> ph) {
    return vprops_.get_property(ph);
  }

  template <class T>
  const property<T>& get_property(vprop_handle<T> ph) const {
    return vprops_.get_property(ph);
  }

  //! TODO: access face/edge/halfedge property

public: // access a property element using a hanlde to a mesh item
  template <class T>
  typename vprop_handle<T>::reference
  get_property(vprop_handle<T> ph, vert_handle vh) {
    return vprops_.get_property(ph)[vh.index()];
  }
  template <class T>
  typename vprop_handle<T>::const_reference
  get_property(vprop_handle<T> ph, vert_handle vh) const {
    return vprops_.get_property(ph)[vh.index()];
  }  

  template <class T>
  typename fprop_handle<T>::reference
  get_property(fprop_handle<T> ph, face_handle fh) {
    return fprops_.get_property(ph)[fh.index()];
  }
  template <class T>
  typename fprop_handle<T>::const_reference
  get_property(fprop_handle<T> ph, face_handle fh) const {
    return fprops_.get_property(ph)[fh.index()];
  }

  template <class T>
  typename eprop_handle<T>::reference
  get_property(eprop_handle<T> ph, edge_handle eh) {
    return eprops_.get_property(ph)[eh.index()];
  }
  template <class T>
  typename eprop_handle<T>::const_reference
  get_property(eprop_handle<T> ph, edge_handle eh) const {
    return eprops_.get_property(ph)[eh.index()];
  }

  template <class T>
  typename hprop_handle<T>::reference
  get_property(hprop_handle<T> ph, halfedge_handle hh) {
    return hprops_.get_property(ph)[hh.index()];
  }

  template <class T>
  typename hprop_handle<T>::const_reference
  get_property(hprop_handle<T> ph, halfedge_handle hh) const {
    return hprops_.get_property(ph)[hh.index()];
  }
  
public:
  /// property number
  size_t n_vprops() const { return vprops_.size(); }
  size_t n_eprops() const { return eprops_.size(); }
  size_t n_fprops() const { return fprops_.size(); }
  size_t n_hprops() const { return hprops_.size(); }

  /// element number for each property
  size_t n_verts() const { return vprops_.n_elements(); }
  size_t n_edges() const { return eprops_.n_elements(); }
  size_t n_faces() const { return fprops_.n_elements(); }
  size_t n_halfedges() const { return hprops_.n_elements(); }

  /// synchronizez properties
  void vprops_resize(size_t n) { vprops_.resize(n); }
  void vprops_reserve(size_t n) { vprops_.reserve(n); }
  void vprops_clear() {vprops_.clear(); }
  void vprops_swap(size_t i0, size_t i1) { vprops_.swap(i0, i1); }
  
  void fprops_resize(size_t n) { fprops_.resize(n); }
  void fprops_reserve(size_t n) { fprops_.reserve(n); }
  void fprops_clear() { fprops_.clear(); }
  void fprops_swap(size_t i0, size_t i1) { fprops_.swap(i0, i1); }
  
  void eprops_resize(size_t n) { eprops_.resize(n); }
  void eprops_reserve(size_t n) { eprops_.reserve(n); }
  void eprops_clear() { eprops_.clear(); }
  void eprops_swap(size_t i0, size_t i1) { eprops_.swap(i0, i1); }
  
  void hprops_resize(size_t n) { hprops_.resize(n); }
  void hprops_reserve(size_t n) { hprops_.reserve(n); }
  void hprops_clear() { hprops_.clear(); }
  void hprops_swap(size_t i0, size_t i1) { hprops_.swap(i0, i1); }

private:  
  property_container vprops_; // vert properties
  property_container fprops_; // face properties
  property_container eprops_; // edge properties
  property_container hprops_; // halfedge properties
};
  
}

#endif
