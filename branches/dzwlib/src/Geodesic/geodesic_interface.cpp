/* read mesh from file and 
 - if one vertex is specified, for all vertices of the mesh print their distances to this vertex
 - if two vertices are specified, print the shortest path between these vertices 

	Danil Kirsanov, 01/2008 
*/
#include <iostream>
#include <fstream>

#include "geodesic_algorithm_exact.h"
#include "geodesic_interface.h"

int get_geodesic_vertex_to_vertex(const std::vector<double>& points,
                                  const std::vector<unsigned>& faces,
                                  unsigned source_vertex, unsigned des_vertex,
                                  std::vector<double>& path,
                                  std::vector<std::pair<unsigned, unsigned> >& edge_pair)
{
  geodesic::Mesh mesh;
  ///create internal mesh data structure including edges
  mesh.initialize_mesh_data(points, faces);
  ///create exact algorithm for the mesh
  geodesic::GeodesicAlgorithmExact algorithm(&mesh);
  ///create source
  geodesic::SurfacePoint source(&mesh.vertices()[source_vertex]);
  ///in general, there could be multiple sources, but now we have only one
  std::vector<geodesic::SurfacePoint> all_sources(1,source);
  ///create target
  geodesic::SurfacePoint target(&mesh.vertices()[des_vertex]);
  ///geodesic path is a sequence of SurfacePoints
  std::vector<geodesic::SurfacePoint> p_path;
  ///find a single source-target path
  algorithm.geodesic(source, target, p_path, edge_pair);
  path.resize(p_path.size() * 3);
  for(unsigned i = 0; i < p_path.size(); ++i) {
    geodesic::SurfacePoint& s = p_path[i];
    path[i * 3] = s.x(); path[i * 3 + 1] = s.y(); path[i * 3 + 2] = s.z();
  }
  return 0;
}
