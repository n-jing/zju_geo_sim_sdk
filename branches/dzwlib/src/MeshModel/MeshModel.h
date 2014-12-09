/**   
 * @file MeshModel.h
 * @brief Define mesh model
 * @author dangzw
 * @date Oct 18, 2011 1:55:47 PM
 * @version V1.0
 */

#ifndef MESHMODEL_H_
#define MESHMODEL_H_

#include <string>

#include <Common/BaseData.h>

namespace dzw {
namespace model {
/** MeshModel class is data structure of the model.
 *  Now the class is simple, include vertex, edge, face, normal, and their relations.
 *  And there are some simple functions about the model.
 *  TODO: add more operations into the class.
 */
class mesh_model {

public:
  mesh_model(const boost::property_tree::ptree& opt_);
  /** Get the vertex coordinate of model.
  * @return[out] vertex coordinate.
  */
  const common::VertexCoordMat& get_vertex_coord_mat() const {
    return vertex_coord;
  }
  /** Get vertex index of model edge.
  * @return[out] vertex index of edge.
  */
  const common::VertexOfEdgeVec& get_edge_vec() const {
    return vertex_of_edge;
  }
  /** Get vertex index of model face.
  * @return[out] vertex index of face.
  */
  const common::VertexOfFaceMat& get_face_mat() const {
    return vertex_of_face;
  }
  /**Get edge index of face.
  * @return[out] edge index of face.
  */
  const common::EdgeOfFaceVec& get_edge_of_face() const {
    return edge_of_face;
  }
  /** Get vertex normal.
  * @return[out] every vertex normal
  */
  const common::VertexNormalMat& get_vertex_normal() const {
    return vertex_normal;
  }
  /** Get face normal.
  * @return[out] every face normal
  */
  const common::FaceNormalMat& get_face_normal() const {
    return face_normal;
  }
  /** Get adjacent face index of vertex.
  * @return[out] adjacent face index of vertex.
  */
  const common::AdjFaceOfVertex& get_adj_face_of_vertex() const {
    return adj_face_of_vertex;
  }
  /** Get adjacent face index of edge.
  * @return[out] adjacent face index of edge.
  */
  const common::AdjFaceOfEdge& get_adj_face_of_edge() const {
    return adj_face_of_edge;
  }
  /** Get adjacent edge index of vertex.
  * @return[out] adjacent edge index of vertex.
  */
  const common::AdjEdgeOfVertex& get_adj_edge_of_vertex() const {
    return adj_edge_of_vertex;
  }

  /** Get file name of model.
  * @return[out] file name, type is string.
  */
  const std::string& get_model_name() const {
    return obj_file_name;
  }
  const std::vector<std::vector<int> >& get_feature_line() const {
    return feature_line_vec;
  }
  /**Get edge index by two vertex index.
  * @param index1 vertex 1.
  * @param index2 vertex 2.
  * @return edge index.
  */
  int get_edge_index(int index1, int index2) const ;
  /**Get opposite face of edge.
  * @param edge_index edge index.
  * @param know_face face index.
  * @return opposite face index, if the edge is boundary return -1.
  */
  int get_opposite_face_of_edge(int edge_index, int know_face) const ;
  /**Get facing vertex of a given edge in a specified face.
   * @param edge_index the given edge index.
   * @param face_index the specified face index.
   * @return the facing vertex index.
   */
  int get_facing_vertex_of_edge(int edge_index, int face_index) const ;
  /**Get facing edge of a given vertex in a specified face.
   * @param vertex_index the given vertex index.
   * @param face_index the specified face index.
   * @return the facing edge index.
   */
  int get_facing_edge_of_vertex(int vertex_index, int face_index) const;
  /**If two edges are in same face.
  * @param index1 edge 1.
  * @param index2 edge 2
  * @return if with same face return the face index, otherwise return -1.
  */
  int face_index_of_two_edge(int index1, int index2) const ;
  /**Get the edge index of two faces.
  * @param index1 face 1.
  * @param index2 face 2.
  * @return if have a common edge return edge index, or return -1.
  */
  int edge_index_of_two_face(int index1, int index2) const ;
  /**Get face area.
  * @return face area.
  */
  const std::vector<double>& get_face_area() const{
    return face_area;
  }
  /**Get bounding sphere info
  * @param center_ bounding sphere center.
  * @param radius_ bounding sphere radius.
  */
  void get_bounding_sphere(common::PointCoord& center_, double& radius_) const{
    center_ = model_bounding_center;
    radius_ = model_bounding_radius;
  }
  /**Get one ring vertex index of one vertex.
  * @param vertex_index_ vertex index.
  * @return one ring vertex index of vertex_index_
  */
  const std::vector<int>& get_one_ring_vertex(int vertex_index_) const {
    return one_ring_vertex_index[vertex_index_];
  }
  /**Get vertex status
   * @return vertex status
   */
  const std::vector<bool>& get_vertex_status() const {
    return vertex_status;
  }
  /**Get edge status
   * @return edge status
   */
  const std::vector<bool>& get_edge_status() const {
    return edge_status;
  }
  /**Get face status
   * @return face status
   */
  const std::vector<bool>& get_face_status() const {
    return face_status;
  }
  /**Whether vertex is boundary.
  * @param index_ vertex index.
  * @return
  */
  bool is_boundary_vertex(int index_) const {
    assert(index_ < vertex_coord.size());
    return boundary_vertex[index_];
  }
  /**Whether edge is boundary.
  * @param index_ edge index.
  * @return
  */
  bool is_boundary_edge(int index_) const {
    assert(index_ < vertex_of_edge.size());
    return boundary_edge[index_];
  }
  /**Whether face is boundary.
  * @param index_ face index.
  * @return
  */
  bool is_boundary_face(int index_) const {
    assert(index_ < vertex_of_face.size());
    return boundary_face[index_];
  }
  /**Whether feature edge
   * @param index_ edge index.
   * @return
   */
  bool is_feature_edge(int v1, int v2) const {
    assert(v1 < vertex_coord.size() && v2 < vertex_coord.size());
    int edge_index = get_edge_index(v1, v2);
    assert(edge_index != -1);
    return feature_edge[edge_index];
  }
  bool is_feature_vertex(int index) const {
    assert(index < vertex_coord.size());
    return feature_vertex[index];
  }
  /**Get dihedral angle of two adjacent faces.
   * @param dihedral_angle returned dihedral angle for each edge.
   * @return
   */
  int get_dihedral_angle(std::vector<double>& dihedral_angle) const ;
  /**Add a vertex on a specified edge.
   * @param edge_index the specified edge index.
   * @param new_vertex the added new vertex coordinate.
   * @return success or fail
   */
  int add_vertex_on_edge(int edge_index, const zjucad::matrix::matrix<double>& new_vertex);
  /**Add a vertex on a specified face.
   * @param face_index the specified face index.
   * @param new_vertex the added new vertex coordinate.
   * @return success or fail
   */
  int add_vertex_on_face(int face_index, const zjucad::matrix::matrix<double>& new_vertex);
  /**Merge a vertex to an edge, the vertex and edge must in same face.
   * @param vertex_index the vertex index
   * @param edge_index the edge index.
   * @param face_index the face index.
   * @return
   */
  int merge_vertex_to_edge(int vertex_index, int edge_index, int face_index);
  /**Merge two close vertex, which must the endpoint of same edge.
   * @param v1 first vertex index.
   * @param v2 second vertex index.
   * @return
   */
  int merge_two_close_vertex(int v1, int v2);
  /**Update a specified vertex coordinate with a given new coordinate.
   * @param index the specified vertex index
   * @param new_coord the new vertex coordinate
   * @return success or fail
   */
  int update_vertex_coord(int index, const zjucad::matrix::matrix<double>& new_coord);
  /**Save mesh model as obj format.
   * @param file_name obj file name.
   * @return success or fail
   */
  int save_mesh_model_as_obj() const ;
  /**Set one feature line.
   * @param one_feature_line one feature line vector.
   */
  void set_one_feature_line(const std::vector<int>& one_feature_line) ;
  /**Save feature lines of mesh model.
    */
  void save_feature_lines();
  /**Set edge(v1, v2) as feature edge.
   * @param v1 one vertex index of feature edge.
   * @param v2 another vertex index of feature edge.
   */
  void set_feature_edge(int v1, int v2) ;
private:
  /**
  * initial mesh model, in this function need to compute some basic information.
  * @param[in] file_name_ file name of the mesh model.
  * @return type is integer, 0 is success, positive integer means different error.
  * @see enum RETURN_TYPE
  */
  int init_model(const std::string& file_name_);
  /** Load obj model, and compute the bounding info of model. Now we only support triangle mesh as input.
  * @param[in] file_name file name of the mesh model.
  * @return type is integer, 0 is success, positive integer means different error.
  * @see enum RETURN_TYPE
  */
  int load_obj_model(const std::string& file_name_);
  /**Load feature line of model
   * @param file_name_ feature line file name.
   * @return
   */
  int load_feature_line(const std::string& file_name_);
  /** Compute face normal.
  */
  void compute_face_normal();
  /** Compute face area.
  */
  void compute_face_area();
  /** Compute vertex normal.
  */
  void compute_vertex_normal();
  /** Find the adjacent info of vertex and edge.
  */
  void find_adj_info();
  /** Find the edge index by two vertex. If new edge it could be added to the end.
  * @param index1 vertex 1
  * @param index2 vertex 2
  * @return return the edge index.
  */
  int find_edge_index(int index1, int index2);
  /**set boundary edge flag, vertex flag, and face flag.
  */
  void find_boundary();
  /**Add a face with three given vertex index.
   * @param face_vertex vertex index of new face.
   * @return success or fail
   */
  int add_face(const std::vector<int>& face_vertex);
  /**Add a face with a given edge index and a given vertex index.
   * @param edge_index the specified edge index.
   * @param vertex_index the specified vertex index.
   * @return success or fail
   */
  int add_face(int edge_index, int vertex_index);
  /**Delete a specified vertex.
   * @param vertex_index the specified vertex index.
   * @return
   */
  int delete_vertex(int vertex_index);
  /**Delete a specified face.
   * @param face_index the specified face index.
   * @return success or fail
   */
  int delete_face(int face_index);
  /**Delete a specified edge.
   * @param edge_index the specified edge index.
   * @return success or fail
   */
  int delete_edge(int edge_index);
private:
  common::VertexCoordMat 		vertex_coord; /**< vertex coordinate. */
  common::VertexOfFaceMat 		vertex_of_face; /**< vertex index of face. */
  common::VertexOfEdgeVec 		vertex_of_edge; /**< vertex index of edge. */
  common::EdgeOfFaceVec			edge_of_face; /**< edge index of face. */
  common::VertexNormalMat 		vertex_normal; /**< normal of vertex. */
  common::FaceNormalMat 		face_normal; /**< normal of face. */

  std::vector<bool>                     vertex_status; /**< mark the vertex whether be removed. */
  std::vector<bool>                     edge_status; /**< mark the edge whether be removed. */
  std::vector<bool>                     face_status; /**< mark the face whether be removed. */

  common::AdjFaceOfVertex 		adj_face_of_vertex; /**< adjacent face index of vertex. */
  common::AdjFaceOfEdge 		adj_face_of_edge; /**< adjacent face index of edge. */
  common::AdjEdgeOfVertex 		adj_edge_of_vertex; /**< adjacent edge index of vertex. */
  std::vector<std::vector<int> > 	one_ring_vertex_index; /**< get one ring vertex index of one vertex. */

  common::PointCoord 			model_bounding_center; /**< the center of model's bounding sphere. */
  double 				model_bounding_radius; /**< the radius of model's bounding sphere. */

  std::vector<double>			face_area; /**< area of face. */
  std::vector<bool>			boundary_edge; /**< whether a boundary edge. */
  std::vector<bool>			boundary_face; /**< whether a boundary face. */
  std::vector<bool>			boundary_vertex; /**< whether a boundary vertex. */
  std::vector<bool>                     feature_edge; /**< whether a feature edge. */
  std::vector<bool>                     feature_vertex; /**< whether a feature vertex. */
  std::vector<std::vector<int> >        feature_line_vec; /**< feature line vector. */

  std::string 				obj_file_name; /**< mesh model file name. */
  std::string                           fl_file_name; /**< feature line file name. */
  std::string                           out_obj_file_name; /**< output mesh model file name. */
  std::string                           out_fl_file_name; /**< output feature line file name. */
  boost::property_tree::ptree 	        opts;
};
}
}
#endif /* MESHMODEL_H_ */
