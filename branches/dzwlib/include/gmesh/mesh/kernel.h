#ifndef GMESH_MESH_KERNEL_H_
#define GMESH_MESH_KERNEL_H_

#include "mesh_item.h"
#include "property_mananger.h"
#include "geometry.h"
#include "halfedge_connectivity.h"

#include <iostream>

#include <vector>
#include <string>
#include <cstdio>

namespace gmesh {

  template<class geometry_t=geometry, class connectivity_t=halfedge_connectivity>
      class kernel : public geometry_t, public connectivity_t {
 public:
 typedef typename connectivity_t::vert_iter vert_iter;
 typedef typename connectivity_t::edge_iter edge_iter;
 typedef typename connectivity_t::face_iter face_iter;
 typedef typename connectivity_t::const_vert_iter const_vert_iter;
 typedef typename connectivity_t::const_edge_iter const_edge_iter;
 typedef typename connectivity_t::const_face_iter const_face_iter;
  
 typedef typename connectivity_t::vv_iter vv_iter;
 typedef typename connectivity_t::const_vv_iter const_vv_iter;

 private:
 using connectivity_t::new_vert;
 using geometry_t::set_point;
 using geometry_t::get_point;
 using connectivity_t::set_halfedge_handle;
 using connectivity_t::set_vert_handle;
 using connectivity_t::set_edge_handle;
 using connectivity_t::set_face_handle;
 using connectivity_t::set_prev_handle;
 using connectivity_t::set_next_handle;
 using connectivity_t::set_oppo_handle;
 using connectivity_t::new_face;
 using connectivity_t::add_edge;
 using connectivity_t::set_boundary;
 
 public:
 using connectivity_t::add_face;
 using connectivity_t::get_halfedge_handle;
 using connectivity_t::get_to_vert_handle;
 using connectivity_t::get_prev_halfedge_handle;
 using connectivity_t::get_next_halfedge_handle;
 using connectivity_t::get_face_handle;
 using connectivity_t::get_adj_face;
 using connectivity_t::get_adj_vert;
 using connectivity_t::get_adj_edge;
 using connectivity_t::get_verts;
 using connectivity_t::get_edges;
 using connectivity_t::get_vert_num;
 using connectivity_t::get_edge_num;
 using connectivity_t::get_face_num;
 using connectivity_t::delete_face;
 using connectivity_t::delete_vert;
 using connectivity_t::delete_edge;
 using connectivity_t::garbage_collection;
 using connectivity_t::is_flip_ok;
 using connectivity_t::flip_edge;
 using connectivity_t::is_collapse_ok;
 using connectivity_t::collapse_edge;
 using connectivity_t::is_boundary;
 using connectivity_t::is_valid_handle;
 
 /* using connectivity_t::verts_begin; */
 /* using connectivity_t::verts_end; */
 /* using connectivity_t::edges_begin; */
 /* using connectivity_t::faces_begin; */
 /* using connectivity_t::faces_end; */
 /* using connectivity_t::vv_circulator; */
 /* using connectivity_t::cvv_circulator; */

 public:
 kernel();
 virtual ~kernel();
 
 vert_handle add_vert(point p);
 inline void set_coord(vert_handle vh, const point &p)
 {
   set_point(vh, p);
 };
  
 const property_mananger& get_property_mananger() const;
 property_mananger& get_property_mananger();

 point get_coord(vert_handle vh) const;

 edge_handle split_edge(edge_handle eh);

 vert_iter verts_begin() { return connectivity_t::verts_begin(); }
 vert_iter verts_end()   { return connectivity_t::verts_end();   }
 edge_iter edges_begin() { return connectivity_t::edges_begin(); }
 edge_iter edges_end()   { return connectivity_t::edges_end();   }
 face_iter faces_begin() { return connectivity_t::faces_begin(); }
 face_iter faces_end()   { return connectivity_t::faces_end();   }

 const_vert_iter verts_begin() const { return connectivity_t::verts_begin(); }
 const_vert_iter verts_end()   const { return connectivity_t::verts_end();   }
 const_edge_iter edges_begin() const { return connectivity_t::edges_begin(); }
 const_edge_iter edges_end()   const { return connectivity_t::edges_end();   }
 const_face_iter faces_begin() const { return connectivity_t::faces_begin(); }
 const_face_iter faces_end()   const { return connectivity_t::faces_end();   }

 vv_iter vv_circulator(vert_handle vh) {
   return connectivity_t::vv_circulator(vh);
 }
 const_vv_iter cvv_circulator(vert_handle vh) const {
   return connectivity_t::cvv_circulator(vh);
 }

  template<typename T>
      void add_property(vprop_handle<T> &vph, const std::string &name)
  {
    pm_.add_property(vph, name);
  }

  template<typename T>
      void add_property(fprop_handle<T> &fph, const std::string &name)
  {
    pm_.add_property(fph, name);
  }

  template<typename T>
      void add_property(eprop_handle<T> &eph, const std::string &name)
  {
    pm_.add_property(eph, name);
  }

  template<typename T>
      void remove_property(vprop_handle<T> &vph)
  {
    pm_.remove_property(vph);
  }

  template<typename T>
      void remove_property(fprop_handle<T> &fph)
  {
    pm_.remove_property(fph);
  }

  template<typename T>
      void remove_property(eprop_handle<T> &eph)
  {
    pm_.remove_property(eph);
  }

  template<typename T>
      bool get_property_handle(vprop_handle<T> &vph, const std::string &name) const
  {
    return pm_.get_property_handle(vph, name);
  }

  template<typename T>
      bool get_property_handle(fprop_handle<T> &fph, const std::string &name) const
  {
    return pm_.get_property_handle(fph, name);
  }
    
  template<typename T>
      bool get_property_handle(eprop_handle<T> &eph, const std::string &name) const
  {
    return pm_.get_property_handle(eph, name);
  }

  template<typename T>
      typename vprop_handle<T>::reference
      get_property(vprop_handle<T> vph, vert_handle vh)
  {
    assert(is_valid_handle(vh));
    return pm_.get_property(vph, vh);
  }

  template<typename T>
      typename fprop_handle<T>::reference
      get_property(fprop_handle<T> fph, face_handle fh)
  {
    assert(is_valid_handle(fh));
    return pm_.get_property(fph, fh);
  }

  template<typename T>
      typename eprop_handle<T>::reference
      get_property(eprop_handle<T> eph, edge_handle eh)
  {
    assert(is_valid_handle(eh));
    return pm_.get_property(eph, eh);
  }
  
  template<typename T>
      typename vprop_handle<T>::const_reference
      get_property(vprop_handle<T> vph, vert_handle vh) const
  {
    assert(is_valid_handle(vh));
    return pm_.get_property(vph, vh);
  }

  template<typename T>
      typename fprop_handle<T>::const_reference
      get_property(fprop_handle<T> fph, face_handle fh) const
  {
    assert(is_valid_handle(fh));
    return pm_.get_property(fph, fh);
  }

  template<typename T>
      typename eprop_handle<T>::const_reference
      get_property(eprop_handle<T> eph, edge_handle eh) const
  {
    assert(is_valid_handle(eh));
    return pm_.get_property(eph, eh);
  }
 

 private:
 property_mananger pm_;
  
  };

  template<class geometry_t, class connectivity_t>
      kernel<geometry_t, connectivity_t>::kernel() : geometry_t(pm_), connectivity_t(pm_) {
    /// default properties
    /// vertex coordinate property
    pm_.add_property(geometry_t::get_vert_coord_property_handle(), "v:points");

    /// vertex status property
    pm_.add_property(connectivity_t::get_vert_status_property_handle(), "v:status");
    /// face status property
    pm_.add_property(connectivity_t::get_face_status_property_handle(), "f:status");
    /// edge status property
    pm_.add_property(connectivity_t::get_edge_status_property_handle(), "e:status");
    
    geometry_t::set_geometry_base_property();
    connectivity_t::set_conn_base_property();
  }

  template<class geometry_t, class connectivity_t>
      kernel<geometry_t, connectivity_t>::~kernel() {}
  
  template<class geometry_t, class connectivity_t>
      const property_mananger& kernel<geometry_t, connectivity_t>::get_property_mananger() const{
    return pm_;
  }

  template<class geometry_t, class connectivity_t>
      property_mananger& kernel<geometry_t, connectivity_t>::get_property_mananger(){
    return pm_;
  }

  template<class geometry_t, class connectivity_t>
      vert_handle kernel<geometry_t, connectivity_t>::add_vert(point p){
    vert_handle vh = new_vert();
    set_point(vh, p);
    return vh;
  }

  template<class geometry_t, class connectivity_t>
      point kernel<geometry_t, connectivity_t>::get_coord(vert_handle vh) const {
    return get_point(vh);
  }

  template<class geometry_t, class connectivity_t>
      edge_handle kernel<geometry_t, connectivity_t>::split_edge(edge_handle eh)
  {
    assert(is_valid_handle(eh));
    halfedge_handle heh0 = get_halfedge_handle(eh, 0);
    halfedge_handle heh1 = get_halfedge_handle(eh, 1);
    face_handle f0 = get_face_handle(heh0);
    face_handle f1 = get_face_handle(heh1);
    std::vector<edge_handle> vv0, vv1;
    if (is_valid_handle(f0)) vv0 = get_edges(f0);
    if (is_valid_handle(f1)) vv1 = get_edges(f1);

    vert_handle v0 = get_to_vert_handle(heh0);
    vert_handle v1 = get_to_vert_handle(heh1);
    point p0 = get_coord(v0);
    point p1 = get_coord(v1);
    vert_handle va = add_vert((p0+p1)*0.5);

    halfedge_handle hep0 = get_prev_halfedge_handle(heh0);
    halfedge_handle hen1 = get_next_halfedge_handle(heh1);
    halfedge_handle hen0 = get_next_halfedge_handle(heh0);
    halfedge_handle hep1 = get_prev_halfedge_handle(heh1);

    vert_handle v2 = get_to_vert_handle(hen0);
    vert_handle v3 = get_to_vert_handle(hen1);

    // split the edge : eh,  add a vertex va.
    edge_handle ehm = add_edge(v1, va);
    halfedge_handle hehm0 = get_halfedge_handle(ehm, 0);
    halfedge_handle hehm1 = get_halfedge_handle(ehm, 1);

    set_halfedge_handle(va, heh0);
    set_vert_handle(heh1, va);
    set_boundary(hep0);
    set_boundary(hen1);

    set_next_handle(hehm0, heh0);
    set_prev_handle(heh0, hehm0);
    set_next_handle(heh1, hehm1);
    set_prev_handle(hehm1, heh1);

    set_next_handle(hehm1, hen1);
    set_prev_handle(hen1, hehm1);
    set_next_handle(hep0, hehm0);
    set_prev_handle(hehm0, hep0);
    set_vert_handle(hehm0, va);
    set_vert_handle(hehm1, v1);
    set_halfedge_handle(v1, hehm0);
    //adjust_vert_halfedge_handle(v1);
    set_boundary(hehm0);
    set_boundary(hehm1);
  
    if (is_valid_handle(f0)) {
      assert(vv0.size() == 3);
      set_halfedge_handle(f0, hen0);

      edge_handle eha = add_edge(v2, va);
      halfedge_handle heha0 = get_halfedge_handle(eha, 0);
      halfedge_handle heha1 = get_halfedge_handle(eha, 1);
      set_face_handle(heha0, f0);
      //set_boundary(heha1);
      // update the link relation : eha, eh, ehm
      set_next_handle(heha0, heh0);
      set_prev_handle(heh0, heha0);
      set_next_handle(hen0, heha0);
      set_prev_handle(heha0, hen0);
      set_next_handle(heha1, hep0);
      set_prev_handle(hep0, heha1);
      set_next_handle(hehm0, heha1);
      set_prev_handle(heha1, hehm0);
      set_next_handle(hep0, hehm0);
      set_prev_handle(hehm0, hep0);

      // add a new face, the incident edges are hep0, heha1, hehm0
      face_handle fa = new_face();
      set_halfedge_handle(fa, hep0);
      set_face_handle(heha1, fa);
      set_face_handle(hep0, fa);
      set_face_handle(hehm0, fa);
    }

    if (is_valid_handle(f1)) {
      assert(vv1.size() == 3);
      set_halfedge_handle(f1, hep1);
      edge_handle ehb = add_edge(v3, va);
      halfedge_handle hehb0 = get_halfedge_handle(ehb, 0);
      halfedge_handle hehb1 = get_halfedge_handle(ehb, 1);
      set_face_handle(hehb1, f1);

      // update the link relation : ehb, eh, ehm
      set_next_handle(hehb0, hehm1);
      set_prev_handle(hehm1, hehb0);
      set_next_handle(hehm1, hen1);
      set_prev_handle(hen1, hehm1);
      set_next_handle(heh1, hehb1);
      set_prev_handle(hehb1, heh1);
      set_next_handle(hehb1, hep1);
      set_prev_handle(hep1, hehb1);
      set_next_handle(hen1, hehb0);
      set_prev_handle(hehb0, hen1);

      // add a new face, the incident edges are hen1, hehm1, hehb0
      face_handle fb = new_face();
      set_halfedge_handle(fb, hen1);
      set_face_handle(hehm1, fb);
      set_face_handle(hen1, fb);
      set_face_handle(hehb0, fb);
    }

    return ehm;
  }
 
  

} //namespace gmesh

#endif 
