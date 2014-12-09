#ifndef GMESH_MESH_GEOMETRY_H_
#define GMESH_MESH_GEOMETRY_H_

#include "handle.h"
#include "property_handle.h"
#include "property_mananger.h"
#include "mesh_item.h"
#include "../common/types.h"

namespace gmesh {

class geometry {
public:
  geometry(property_mananger& pm) : pm_(pm) {}
  virtual ~geometry() {}

  const point* get_points() const;
  const point& get_point(vert_handle vh) const;
  void set_point(vert_handle vh, const point& p);

  const vprop_handle<point>& get_vert_coord_property_handle() const;
  vprop_handle<point>& get_vert_coord_property_handle();

protected:
  void set_geometry_base_property();
  
private:
  property_mananger& pm_;
  vprop_handle<point> vc_handle_; //vert point property handle
};

inline const point* geometry::get_points() const {
  return pm_.get_property(vc_handle_).data();
}

inline const point& geometry::get_point(vert_handle vh) const {
  return pm_.get_property(vc_handle_, vh);
}

inline void geometry::set_point(vert_handle vh, const point& p){
  pm_.get_property(vc_handle_, vh) = p;
}

inline const vprop_handle<point>& geometry::get_vert_coord_property_handle() const{
  return vc_handle_;
}

inline vprop_handle<point>& geometry::get_vert_coord_property_handle(){
  return vc_handle_;
}

inline void geometry::set_geometry_base_property()
{
  
}

}

#endif

