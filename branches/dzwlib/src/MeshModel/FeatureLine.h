/**   
* @file FeatureLine.h
* @brief Define feature line
* @author dangzw
* @date Oct 20, 2011 4:56:41 PM 
* @version V1.0   
*/

#ifndef FEATURELINE_H_
#define FEATURELINE_H_

#include <boost/shared_ptr.hpp>
#include <vector>
#include <Common/Complex.h>
#include "MeshModel.h"

namespace dzw{
namespace model{

/**Feature line of the model, which can be used to express the hard constrain.
 * We also can convert the constrain edge into constrain face.
 */
class feature_line{
public:
  /**Constructor
  * @return
  */
  feature_line(){}
  /**Constructor
  * @param mesh_model_
  * @return
  */
  feature_line(const boost::property_tree::ptree& opt_);
  /**Destructor.
  * @return
  */
  ~feature_line(){}
  /**Get the feature line vertex coordinate.
   * @return the vertex coordinate.
   */
  const dzw::common::VertexCoordMat& get_feature_point() const {
    return feature_point_coord;
  }
  /**Get the feature line edges
   * @return the vector of feature edges.
   */
  const std::vector<std::vector<int> >& get_feature_edge() const {
    return feature_edge_vec;
  }
  /**Get edge index with given vertex 1 and vertex 2.
   * @param vertex_1 vertex 1 index
   * @param vertex_2 vertex 2 index
   * @return edge index
   */
  int get_edge_index(int vertex_1, int vertex_2) const ;
private:
  /**Load feature line from *.fl file.
  * @param file_name feature line file name.
  * @return 0 is success, others stand for different error.
  */
  int load_feature_line(const std::string& file_name);
  /**Compute model bounding box and bounding radius.
   */
  void compute_bounding_box();
  /**Find adjacent information of basic element of mesh.
   */
  void find_adjacent_info();

private:
  dzw::common::VertexCoordMat                           feature_point_coord;
  std::vector<std::vector<int> >                        feature_edge_vec;

  std::vector<std::vector<int> >                        adj_edge_index_of_vertex;

  boost::property_tree::ptree 	                        opts;
  common::PointCoord                                    model_bounding_center; /**< the center of model's bounding sphere. */
  double                                                model_bounding_radius; /**< the radius of model's bounding sphere. */
};
}
}
#endif /* FEATURELINE_H_ */
