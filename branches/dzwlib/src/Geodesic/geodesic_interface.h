/**   
* @file geodesic_interface.h 
* @brief TODO
* @author dangzw
* @date Apr 27, 2012 8:12:45 PM 
* @version V1.0   
*/

#ifndef GEODESIC_INTERFACE_H_
#define GEODESIC_INTERFACE_H_

#include <Common/BaseData.h>

/** Get the geodesic path between two vertex.
  * @param points mesh vertex
  * @param faces face vector
  * @param source_vertex source vertex index
  * @param des_vertex destination vertex index
  * @param path the passing points coordinate
  * @param edge_index the pasing edge index vector
  */

int get_geodesic_vertex_to_vertex(const std::vector<double>& points,
                                  const std::vector<unsigned>& faces,
                                  unsigned source_vertex, unsigned des_vertex,
                                  std::vector<double>& path,
                                  std::vector<std::pair<unsigned, unsigned> >& edge_index);


#endif /* GEODESIC_INTERFACE_H_ */
