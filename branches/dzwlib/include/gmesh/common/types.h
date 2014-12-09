#ifndef MESHLIB_BASICDATATYPE_H_

#include "macro.h"
#include "vec2d.h"
#include "vec3d.h"
#include "color.h"

#include <vector>

namespace gmesh{
  /* ================== Various Indices ================== */
  typedef int VertexID;    // Vertex index
  typedef int	EdgeID;      // HalfEdge index
  typedef int	FaceID;      // Face index

  typedef int	PointID;     // Point index
  typedef int	LineID;      // Line index
  typedef int	PolygonID;   // Polygon index

  typedef int BdyID;       // Boundary index
  typedef int ComponentID; // Component index

  typedef int	ModelID;     // Model index
  typedef int	SceneID;     // Scene index
    
  /* ================== Properties ================== */
  typedef Vec3D<double> Coord;
  typedef Vec3D<double> Coord3D;
  typedef Vec2D<double> Coord2D;
    
  typedef Coord   Normal;
  typedef Coord2D TexCoord;
  typedef int     Flag;
  typedef int     Index;
  
  typedef unsigned int uint;
  typedef Coord point;
  typedef Coord vec3d;
    
} // namespace meshlib

#endif
