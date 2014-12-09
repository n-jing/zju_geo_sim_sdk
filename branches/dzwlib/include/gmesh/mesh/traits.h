#ifndef GMESH_MESH_TRAITS_H_
#define GMESH_MESH_TRAITS_H_

#include "../common/types.h"

namespace gmesh {

struct default_traits {
  typedef vec3f point; // default coordinate type is vec3f
  typedef vec3f normal;
  typedef float texcoord1d;
  typedef vec2f texcoord2d;
  typedef vec3f texcoord3d;
};

}

#endif
