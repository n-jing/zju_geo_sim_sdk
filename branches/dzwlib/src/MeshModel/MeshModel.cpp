

/**
* @file mesh_model.cpp
* @brief TODO
* @author dangzw
* @date Oct 21, 2011 4:03:59 PM
* @version V1.0
*/
#include <fstream>
#include <cassert>
#include <zjucad/matrix/io.h>
#include <boost/lexical_cast.hpp>

#include "../Common/MyException.h"
#include "../Common/Complex.h"
#include "MeshModel.h"

using namespace std;
using namespace zjucad::matrix;
using namespace dzw::common;

dzw::model::mesh_model::mesh_model(const boost::property_tree::ptree& opt_) : opts(opt_){
  obj_file_name = opts.get<string>("input/model_file_name.value");
  int prefix_len = obj_file_name.rfind('.');
  std::string suffix = obj_file_name.substr(prefix_len + 1, obj_file_name.length() - prefix_len - 1);
  if(zjucad::has("output/model_file_name.value", opts))
    out_obj_file_name = opts.get<string>("output/model_file_name.value");
  else
    out_obj_file_name = suffix + "_out.obj";
  bool with_fl = opts.get<string>("flag/feature_line.value", "no") == "yes" ? true : false;
  if(with_fl && zjucad::has("input/feature_line.value", opts))
    fl_file_name = opts.get<string>("input/feature_line.value");
  else if(with_fl && !zjucad::has("input/feature_line.value", opts))
    fl_file_name = suffix + ".fl";
  else
    fl_file_name = "";
  if(with_fl && zjucad::has("output/feature_line.value", opts))
    out_fl_file_name = opts.get<string>("output/feature_line.value");
  else if(with_fl && !zjucad::has("output/feature_line.value", opts))
    out_fl_file_name = suffix + "_out.fl";
  init_model(obj_file_name);
  if(with_fl)
    load_feature_line(fl_file_name);
}

int dzw::model::mesh_model::init_model(const std::string& file_name_){
  int prefix_len = obj_file_name.rfind('.');
  std::string suffix = obj_file_name.substr(prefix_len + 1, obj_file_name.length() - prefix_len - 1);
#ifdef _DEBUG
  cout << "[Debug Info] : " << "file name is " << file_name << endl;
  cout << "[Debug Info] : " << "prefix length is " << prefix_len << endl;
  cout << "[Debug Info] : " << "model file format is " << suffix << endl;
#endif
  if(suffix == "obj")
    load_obj_model(file_name_);
  else
    throw my_exception("No support " + suffix + "format.");
#ifdef _DEBUG
  cout << "[Debug info] : load obj file success!" << endl;
#endif
  find_adj_info();
#ifdef _DEBUG
  cout << "[Debug info] : find adjacent info success!" << endl;
#endif
  find_boundary();
#ifdef _DEBUG
  cout << "[Debug info] : find boundary info success!" << endl;
#endif
  compute_face_normal();
#ifdef _DEBUG
  cout << "[Debug info] : compute face normal success!" << endl;
#endif
  compute_vertex_normal();
#ifdef _DEBUG
  cout << "[Debug info] : compute vertex normal success!" << endl;
#endif
//  order_adj_edge_of_vertex();
#ifdef _DEBUG
  cout << "[Debug info] : order adjacent edge of vertex success!" << endl;
#endif
//  order_adj_face_of_vertex();
#ifdef _DEBUG
  cout << "[Debug info] : order adjacent face of vertex success!" << endl;
#endif
  compute_face_area();
#ifdef _DEBUG
  cout << "[Debug info] : compute face area success!" << endl;
#endif
//  compute_angle_defect();
#ifdef _DEBUG
  cout << "[Debug info] : compute angle defect success!" << endl;
#endif
  return 0;
}

int dzw::model::mesh_model::load_obj_model(const std::string & file_name_) {
#ifdef _DEBUG
  cout << "[Debug Info] : obj file name is " << file_name_ << endl;
#endif
  ifstream ifs(file_name_.c_str());
  if(!ifs)
    throw my_exception("Can not open the obj model.");
  int nVertex = 0, nFace = 0; /// Vertex number and face number.
  matrix<double> point_pos(3,1);
  char format[10];
  vector<int> face_i;
  std::string str;
  while(getline(ifs, str)) {
    if(str[0] == '#')	continue;
    if(str[0] == 'g')	continue;
    if(str[0] == 'v' && str[1] == ' ') {
      ++nVertex;
      continue;
    }
    if(str[0] == 'f' && str[1] == ' ') {
      ++nFace;
      continue;
    }
  }
  ifs.clear();
  ifs.seekg(0, std::ifstream::beg);

  /// Store vertex information
  vertex_coord.resize(nVertex);
  vertex_status.resize(nVertex, true);
  /// Store face information
  vertex_of_face.resize(nFace);
  face_status.resize(nFace, true);
  size_t vn = 0, fn = 0;
  while(getline(ifs, str)) {
    if(str[0] == '#')	continue;
    if(str[0] == 'g')	continue;
    if(str[0] == 'o')	continue;
    // read the vertex information here.
    if(str[0] == 'v' && str[1] == ' ') {
      vertex_coord[vn].resize(3, 1);
      sscanf(str.c_str(), "%s %lf %lf %lf", format, &vertex_coord[vn][0], &vertex_coord[vn][1], &vertex_coord[vn][2]);
      ++vn;
      continue;
    }
    // read the face information here.
    if(str[0] == 'f' && str[1] == ' ') {
      std::vector<int> word_pos(1, 0);
      std::string delimit = " \t";
      size_t str_len = str.length();
      for(size_t i = 0; i < str_len; ++i) {
        if(delimit.find(str[i]) == std::string::npos)
          continue;
        if(i+1 < str_len && delimit.find(str[i+1]) == std::string::npos)
          word_pos.push_back((int) i+1);
      }
      if(word_pos.size() != 4) {
        throw my_exception("Our input must be a triangle mesh.");
      }
      vertex_of_face[fn].resize(3, 1);
      sscanf(str.c_str(), "%s %d %d %d", format, &vertex_of_face[fn][0], &vertex_of_face[fn][1], &vertex_of_face[fn][2]);
      vertex_of_face[fn][0]--; vertex_of_face[fn][1]--; vertex_of_face[fn][2]--;
      ++fn;
      continue;
    }
  }
#ifdef _DEBUG
  cout << "[Debug Info] : load obj file over!" << endl;
#endif
  model_bounding_center.resize(3, 1);
  matrix<double> max_point = ones<double>(3, 1) * std::numeric_limits<double>::infinity() * -1.0;
  matrix<double> min_point = ones<double>(3, 1) * std::numeric_limits<double>::infinity();
  for(size_t i = 0; i < 3; ++i){
    for(size_t j = 0; j < vertex_coord.size(); ++j) {
      if(max_point(i, 0) < vertex_coord[j][i])
        max_point(i, 0) = vertex_coord[j][i];
      if(min_point(i, 0) > vertex_coord[j][i])
        min_point(i, 0) = vertex_coord[j][i];
    }
  }
  model_bounding_center = (max_point + min_point) / 2;
  model_bounding_radius = norm(max_point - min_point) / 2;
  ifs.close();
  return 0;
}

int dzw::model::mesh_model::load_feature_line(const std::string& file_name_) {
  std::ifstream ifs(file_name_.c_str());
  int fl_num = 0;
  ifs >> fl_num;
  feature_line_vec.resize(fl_num);
  for(size_t i = 0; i < fl_num; ++i) {
    int vertex_num = 0;
    ifs >> vertex_num;
    feature_line_vec[i].resize(vertex_num);
    for(size_t j = 0; j < vertex_num; ++j) {
      ifs >> feature_line_vec[i][j];
    }
  }
  for(size_t i = 0; i < feature_line_vec.size(); ++i) {
    if(feature_line_vec[i].size() == 1)
      feature_vertex[feature_line_vec[i][0]] = true;
    else {
      for(size_t j = 1; j < feature_line_vec[i].size(); ++j) {
        int v1 = feature_line_vec[i][j - 1], v2 = feature_line_vec[i][j];
        int edge_index = get_edge_index(v1, v2);
        assert(edge_index != -1);
        feature_vertex[v1] = feature_vertex[v2] = true;
        feature_edge[edge_index] = true;
      }
    }
  }
  return 0;
}

int dzw::model::mesh_model::save_mesh_model_as_obj() const {
  ofstream ofs(out_obj_file_name.c_str());
  if(!ofs)
    throw my_exception("Can not open the obj model.");
  std::vector<int> real_index(vertex_coord.size(), -1);
  if(vertex_status[0])
    real_index[0] = 0;
  for(size_t i = 1; i < vertex_status.size(); ++i) {
    if(vertex_status[i])
      real_index[i] = real_index[i - 1] + 1;
    else
      real_index[i] = real_index[i - 1];
  }
  for(size_t i = 0; i < vertex_coord.size(); ++i) {
    if(vertex_status[i]) {
      ofs << "v";
      for(size_t j = 0; j < vertex_coord[i].size(1); ++j)
        ofs << " " << vertex_coord[i][j];
      ofs << endl;
    }
  }
  for(size_t i = 0; i < vertex_of_face.size(); ++i) {
    if(face_status[i]) {
      ofs << "f";
      for(size_t j = 0; j < vertex_of_face[i].size(); ++j)
        ofs << " " << real_index[vertex_of_face[i][j]] + 1;
      ofs << endl;
    }
  }
  ofs.close();
  return 0;
}

void dzw::model::mesh_model::save_feature_lines() {
  std::ofstream ofs(out_fl_file_name.c_str());
  ofs << feature_line_vec.size() << std::endl;
  for(size_t i = 0; i < feature_line_vec.size(); ++i) {
    ofs << feature_line_vec[i].size() << std::endl;
    for(size_t j = 0; j < feature_line_vec[i].size(); ++j)
      ofs << feature_line_vec[i][j] << " ";
    ofs << std::endl;
  }
  ofs.close();
}


int dzw::model::mesh_model::get_edge_index(int index1, int index2) const {
  if(index1 >= vertex_coord.size() || index2 >= vertex_coord.size() || index2 < 0 || index1 < 0)
    throw my_exception("Vertex index out of the mesh vertex number.");
  const std::vector<int>& edge1_vec = adj_edge_of_vertex[index1]; /// we get the adjacent edge of vertex 1.
  const std::vector<int>& edge2_vec = adj_edge_of_vertex[index2]; /// we get the adjacent edge of vertex 2.
  /// find the same edge index between two vector, if exist return the index, else return error.
  for(size_t i = 0; i < edge1_vec.size(); ++i){
    for(size_t j = 0; j < edge2_vec.size(); ++j){
      if(edge1_vec[i] == edge2_vec[j] && edge_status[edge1_vec[i]])
        return edge1_vec[i]; /// successfully find the edge index.
    }
  }
  return -1;
}

int dzw::model::mesh_model::find_edge_index(int index1, int index2){
  if(index1 >= vertex_coord.size() || index2 >= vertex_coord.size() || index1 == index2 || index2 < 0 || index1 < 0)
    throw my_exception("Vertex index out of range or two index is same in function find_edge_index.");
  if(!vertex_status[index1] || !vertex_status[index2])
    throw my_exception("Vertex has already been removed!");
  const std::vector<int>& edge1_vec = adj_edge_of_vertex[index1]; /// we get the adjacent edge of vertex 1.
  const std::vector<int>& edge2_vec = adj_edge_of_vertex[index2]; /// we get the adjacent edge of vertex 2.
  /// find the same edge index between two vector, if exist return the index, else add the new edge.
  for(size_t i = 0; i < edge1_vec.size(); ++i){
    for(size_t j = 0; j < edge2_vec.size(); ++j){
      if(edge1_vec[i] == edge2_vec[j] && edge_status[edge1_vec[i]])
        return edge1_vec[i]; /// successfully find the edge index.
    }
  }
  vector<int> new_edge(2); /// new edge.
  new_edge[0] = index1;
  new_edge[1] = index2;
  vertex_of_edge.push_back(new_edge);
  vector<int> new_edge_null;
  adj_face_of_edge.push_back(new_edge_null);
  edge_status.push_back(true);
  feature_edge.push_back(false);
  one_ring_vertex_index[index1].push_back(index2);
  one_ring_vertex_index[index2].push_back(index1);
  return vertex_of_edge.size() - 1;
}

void dzw::model::mesh_model::find_adj_info() {
  edge_of_face.resize(vertex_of_face.size());
  adj_face_of_vertex.resize(vertex_coord.size());
  adj_edge_of_vertex.resize(vertex_coord.size());
  one_ring_vertex_index.resize(vertex_coord.size());
  /// for each face find the adjacent information.
  for(size_t i = 0; i < vertex_of_face.size(); ++i){
    const vector<int>& face_i = vertex_of_face[i];
    /// for each vertex of face i find the adjacent information.
    for(size_t j = 0; j < face_i.size(); ++j){
      if(face_i[j] >= vertex_coord.size())
        throw my_exception("Face vertex index out of range in finding adjacent info.");
      adj_face_of_vertex[face_i[j]].push_back(i); /// add face i into the adjacent face vector of vertex face[i][j].
      /// find the edge between vertex face[i][j] and face[i][(j+1)%3], if not exist, add the edge.
      int new_edge_index = find_edge_index(face_i[j], face_i[(j + 1) % face_i.size()]);
      edge_of_face[i].push_back(new_edge_index); /// add edge index to face vector.
      if(find(adj_edge_of_vertex[face_i[j]].begin(), adj_edge_of_vertex[face_i[j]].end(), new_edge_index) == adj_edge_of_vertex[face_i[j]].end())
        adj_edge_of_vertex[face_i[j]].push_back(new_edge_index); /// add adjacent edge index of vertex face_i[i][j].
      if(find(adj_edge_of_vertex[face_i[(j + 1) % face_i.size()]].begin(), adj_edge_of_vertex[face_i[(j + 1) % face_i.size()]].end(), new_edge_index)==
          adj_edge_of_vertex[face_i[(j + 1) % face_i.size()]].end())
        adj_edge_of_vertex[face_i[(j + 1) % face_i.size()]].push_back(new_edge_index); /// add adjacent edge index of vertex face_i[i][(j + 1) % face_i.size(1)].
      if(adj_face_of_edge.size() <= new_edge_index){
        vector<int> new_adj_face_of_edge;
        adj_face_of_edge.push_back(new_adj_face_of_edge);
      }
      if(new_edge_index >= adj_face_of_edge.size())
        throw("Edge index out of range in finding adjacent info.");
      adj_face_of_edge[new_edge_index].push_back(i); /// add adjacent face index of edge.
    }
  }
}

void dzw::model::mesh_model::compute_face_normal(){
  face_normal.resize(vertex_of_face.size());
#ifdef _DEBUG
  cout << "[Debug info] : face norma size is " << face_normal.size(2) << endl;
#endif
  for(size_t i = 0; i < vertex_of_face.size(); ++i){
    ///Get the triangle coordinates, this is a 3*3 matrix.
    std::vector<matrix<double> > triangle_coord(vertex_of_face[i].size());
    for(size_t j = 0; j < vertex_of_face[i].size(); ++j) {
      triangle_coord[j] = vertex_coord[vertex_of_face[i][j]];
    }
    dzw::common::normal(triangle_coord, face_normal[i]);
  }
}

void dzw::model::mesh_model::compute_vertex_normal(){
  vertex_normal.resize(vertex_coord.size());
  for(size_t i = 0; i < vertex_normal.size(); ++i){
    vertex_normal[i] = zeros<double>(3, 1);
    for(size_t j = 0; j < adj_face_of_vertex[i].size(); ++j){
      vertex_normal[i] += face_normal[adj_face_of_vertex[i][j]] / norm(face_normal[adj_face_of_vertex[i][j]]);
    }
    vertex_normal[i] /= norm(vertex_normal[i]);
  }
}
void dzw::model::mesh_model::compute_face_area(){
  face_area.resize(vertex_of_face.size());
  for(size_t i = 0; i < vertex_of_face.size(); ++i){
    ///Triangle coordinates.
    std::vector<matrix<double> > triangle_coord(vertex_of_face[i].size());
    for(size_t j = 0; j < vertex_of_face[i].size(); ++j) {
      triangle_coord[j] = vertex_coord[vertex_of_face[i][j]];
    }
    face_area[i] = dzw::common::area(triangle_coord);
    if(face_area[i] < nearly_zero * nearly_zero) {
      std::cout << vertex_of_face[i][0] + 1 << " " << vertex_of_face[i][1] + 1 << " " << vertex_of_face[i][2] + 1 << std::endl;
      throw my_exception("Face area is zero.");
    }
  }
}

int dzw::model::mesh_model::get_opposite_face_of_edge(int edge_index, int know_face) const{
  if(edge_index >= adj_face_of_edge.size() || edge_index < 0)
    throw my_exception("Edge index is out of range in getting opposite face.");
  if(!edge_status[edge_index] || !face_status[know_face])
    throw my_exception("Edge or face has already been removed!");
  int _flag = -1;
  for(size_t i = 0; i < adj_face_of_edge[edge_index].size(); ++i){
    if(adj_face_of_edge[edge_index][i] == know_face)
      _flag = (i + 1) % adj_face_of_edge[edge_index].size();
  }
  if(_flag == -1)
    throw my_exception("Given face is not an adjacent face of edge index. ");
  if(adj_face_of_edge[edge_index].size() == 1)
    return -1;
  return adj_face_of_edge[edge_index][_flag];
}

int dzw::model::mesh_model::get_facing_vertex_of_edge(int edge_index, int face_index) const {
  assert(vertex_of_edge[edge_index].size() == 2);
  for(size_t i = 0; i < vertex_of_face[face_index].size(); ++i) {
    int v_index = vertex_of_face[face_index][i];
    if(vertex_of_edge[edge_index][0] != v_index && vertex_of_edge[edge_index][1] != v_index)
      return v_index;
  }
  return -1;
}

int dzw::model::mesh_model::get_facing_edge_of_vertex(int vertex_index, int face_index) const {
  for(size_t i = 0; i < vertex_of_face[face_index].size(); ++i) {
    if(vertex_of_face[face_index][i] == vertex_index) {
      return get_edge_index(vertex_of_face[face_index][(i + 1) % 3], vertex_of_face[face_index][(i + 2) % 3]);
    }
  }
  return -1;
}
void dzw::model::mesh_model::find_boundary(){
  boundary_edge.resize(vertex_of_edge.size(), false);
  boundary_face.resize(vertex_of_face.size(), false);
  boundary_vertex.resize(vertex_coord.size(), false);
  feature_vertex.resize(vertex_coord.size(), false);
  feature_edge.resize(vertex_of_edge.size(), false);
  for(size_t i = 0; i < vertex_of_edge.size(); ++i){
    if(adj_face_of_edge[i].size() == 1){
      boundary_edge[i] = true;
      boundary_face[adj_face_of_edge[i][0]] = true;
      for(size_t j = 0; j < vertex_of_edge[i].size(); ++j){
        boundary_vertex[vertex_of_edge[i][j]] = true;
      }
    }
    else if(adj_face_of_edge[i].size() > 2)
      throw my_exception("In the mesh, there is a bad edge, which has three adjacent faces.");
  }
}

int dzw::model::mesh_model::face_index_of_two_edge(int index1, int index2) const{
  if(index1 >= adj_face_of_edge.size() || index2 >= adj_face_of_edge.size() || index2 < 0 || index1 < 0)
    throw my_exception("Edge index is out of range in getting face index of two edges.");
  if(!edge_status[index1] || !edge_status[index2])
    throw my_exception("Edge has already been removed!");
  const std::vector<int>& adj_face_1 = adj_face_of_edge[index1];
  const std::vector<int>& adj_face_2 = adj_face_of_edge[index2];
  for(size_t i = 0; i < adj_face_1.size(); ++i){
    for(size_t j = 0; j < adj_face_2.size(); ++j){
      if(adj_face_1[i] == adj_face_2[j])
        return adj_face_1[i];
    }
  }
  return -1;
}

int dzw::model::mesh_model::edge_index_of_two_face(int index1, int index2) const{
  if(index1 >= edge_of_face.size() || index2 >= edge_of_face.size() || index1 < 0 || index2 < 0)
    throw my_exception("Face index is out of range in getting edge index of two faces.");
  if(!face_status[index1] || !face_status[index2])
    throw my_exception("Face has already been removed!");
  const std::vector<int>& adj_edge_1 = edge_of_face[index1];
  const std::vector<int>& adj_edge_2 = edge_of_face[index2];
  for(size_t i = 0; i < adj_edge_1.size(); ++i){
    for(size_t j = 0; j < adj_edge_2.size(); ++j){
      if(adj_edge_1[i] == adj_edge_2[j])
        return adj_edge_1[i];
    }
  }
  return -1;
}

int dzw::model::mesh_model::get_dihedral_angle(std::vector<double>& dihedral_angle) const {
  dihedral_angle.resize(adj_face_of_edge.size(), 0.0);
  for(size_t i = 0; i < adj_face_of_edge.size(); ++i) {
    if(adj_face_of_edge[i].size() == 2) {
      int face_i = adj_face_of_edge[i][0], face_j = adj_face_of_edge[i][1];
      dihedral_angle[i] = dzw::common::angle(face_normal[face_i], face_normal[face_j]);
    }
  }
  return 0;
}

int dzw::model::mesh_model::add_vertex_on_edge(int edge_index, const zjucad::matrix::matrix<double>& new_vertex) {
  vertex_coord.push_back(new_vertex);
  vertex_status.push_back(true);
  std::vector<int> new_null_vector;
  zjucad::matrix::matrix<double> new_null_matrix;
  vertex_normal.push_back(new_null_matrix);
  adj_face_of_vertex.push_back(new_null_vector);
  adj_edge_of_vertex.push_back(new_null_vector);
  one_ring_vertex_index.push_back(new_null_vector);
  boundary_vertex.push_back(false);
  feature_vertex.push_back(false);

//  std::vector<int> new_face_index(3);
  std::vector<int> delete_face_list;
  for(size_t i = 0; i < adj_face_of_edge[edge_index].size(); ++i) {
    int face_index = adj_face_of_edge[edge_index][i];
    for(size_t j = 0; j < edge_of_face[face_index].size(); ++j) {
      if(edge_index != edge_of_face[face_index][j]){
        std::vector<int> new_face(3);
        new_face[0] = vertex_of_face[face_index][j];
        new_face[1] = vertex_of_face[face_index][(j + 1) % vertex_of_face[face_index].size()];
        new_face[2] = vertex_coord.size() - 1;
        add_face(new_face);
//        add_face(edge_of_face[face_index][j], vertex_coord.size() - 1);
      }
    }
    delete_face_list.push_back(face_index);
  }
  for(size_t i = 0; i < delete_face_list.size(); ++i)
    delete_face(delete_face_list[i]);
  delete_edge(edge_index);
  return vertex_coord.size() - 1;
}

int dzw::model::mesh_model::add_vertex_on_face(int face_index, const zjucad::matrix::matrix<double>& new_vertex) {
  vertex_coord.push_back(new_vertex);
  vertex_status.push_back(true);
  std::vector<int> new_null_vector;
  zjucad::matrix::matrix<double> new_null_matrix;
  vertex_normal.push_back(new_null_matrix);
  adj_face_of_vertex.push_back(new_null_vector);
  adj_edge_of_vertex.push_back(new_null_vector);
  one_ring_vertex_index.push_back(new_null_vector);
  boundary_vertex.push_back(false);
  feature_vertex.push_back(false);

  std::vector<int> new_face_index(3);
  for(size_t i = 0; i < edge_of_face[face_index].size(); ++i) {
    std::vector<int> new_face(3);
    new_face[0] = vertex_of_face[face_index][i];
    new_face[1] = vertex_of_face[face_index][(i + 1) % vertex_of_face[face_index].size()];
    new_face[2] = vertex_coord.size() - 1;
    add_face(new_face);
//    add_face(edge_of_face[face_index][i], vertex_coord.size() - 1);
  }
  delete_face(face_index);
  return vertex_coord.size() - 1;
}

int dzw::model::mesh_model::merge_vertex_to_edge(int vertex_index, int edge_index, int face_index) {
  int opp_face_index = get_opposite_face_of_edge(edge_index, face_index);
  if(opp_face_index == -1)
    delete_face(face_index);
  else {
    delete_face(face_index);
    for(size_t i = 0; i < edge_of_face[opp_face_index].size(); ++i) {
      if(edge_of_face[opp_face_index][i] != edge_index) {
        std::vector<int> new_face(3);
        new_face[0] = vertex_of_face[opp_face_index][i];
        new_face[1] = vertex_of_face[opp_face_index][(i + 1) % vertex_of_face[opp_face_index].size()];
        new_face[2] = vertex_coord.size() - 1;
        add_face(new_face);
//        add_face(edge_of_face[opp_face_index][i], vertex_index);
      }
    }
    delete_face(opp_face_index);
  }
  return 0;
}

int dzw::model::mesh_model::merge_two_close_vertex(int v1, int v2) {
  int edge_index = get_edge_index(v1, v2);
  assert(edge_index >= 0);
  bool v1_flag = is_feature_vertex(v1);
  bool v2_flag = is_feature_vertex(v2);
  if(v1_flag && !v2_flag) {
    for(size_t i = 0; i < adj_face_of_vertex[v2].size(); ++i) {
      int candidate_edge_index = get_facing_edge_of_vertex(v2, adj_face_of_vertex[v2][i]);
      assert(vertex_of_edge[candidate_edge_index].size() == 2);
      if(vertex_of_edge[candidate_edge_index][0] != v1 && vertex_of_edge[candidate_edge_index][1] != v1) {
        std::vector<int> new_face = vertex_of_face[adj_face_of_vertex[v2][i]];
        for(size_t j = 0; j < new_face.size(); ++j) {
            if(new_face[j] == v2)
              new_face[j] = v1;
        }
        add_face(new_face);
//        add_face(candidate_edge_index, v1);
      }
    }
    delete_vertex(v2);
  }
  else if(!v1_flag) {
    for(size_t i = 0; i < adj_face_of_vertex[v1].size(); ++i) {
      int candidate_edge_index = get_facing_edge_of_vertex(v1, adj_face_of_vertex[v1][i]);
      assert(vertex_of_edge[candidate_edge_index].size() == 2);
      if(vertex_of_edge[candidate_edge_index][0] != v2 && vertex_of_edge[candidate_edge_index][1] != v2) {
        std::vector<int> new_face = vertex_of_face[adj_face_of_vertex[v1][i]];
        for(size_t j = 0; j < new_face.size(); ++j) {
            if(new_face[j] == v1)
              new_face[j] = v2;
        }
//        add_face(candidate_edge_index, v2);
      }
    }
    delete_vertex(v1);
  }
  else {
    assert(0);
  }
  return 0;
}

int dzw::model::mesh_model::add_face(const std::vector<int>& face_vertex) {
///TODO
  std::vector<int> new_face = face_vertex;
  ///add a new face.
  vertex_of_face.push_back(new_face);
  face_status.push_back(true);
  int new_face_index = vertex_of_face.size() - 1;
  ///add new edges.
  std::vector<int> new_edge_of_face(3);
  for(size_t i = 0; i < new_face.size(); ++i) {
    int new_edge_index = find_edge_index(new_face[i], new_face[(i + 1) % new_face.size()]);
    new_edge_of_face[i] = new_edge_index;
  }
  ///add edge index for new face.
  edge_of_face.push_back(new_edge_of_face);
  ///compute face normal and face area.
  std::vector<matrix<double> > triangle_coord(vertex_of_face[new_face_index].size());
  for(size_t j = 0; j < triangle_coord.size(); ++j) {
    triangle_coord[j] = vertex_coord[vertex_of_face[new_face_index][j]];
  }
  zjucad::matrix::matrix<double> new_face_normal;
  dzw::common::normal(triangle_coord, new_face_normal);
  face_normal.push_back(new_face_normal);
  double new_face_area = dzw::common::area(triangle_coord);
  face_area.push_back(new_face_area);
  ///update adjacent information.
  for(size_t i = 0; i < new_face.size(); ++i) {
    adj_face_of_vertex[new_face[i]].push_back(new_face_index);
  }
  for(size_t i = 0; i < new_face.size(); ++i) {
    if(find(adj_edge_of_vertex[new_face[i]].begin(), adj_edge_of_vertex[new_face[i]].end(), new_edge_of_face[(i + 2) % new_face.size()]) ==
        adj_edge_of_vertex[new_face[i]].end()) {
      adj_edge_of_vertex[new_face[i]].push_back(new_edge_of_face[(i + 2) % new_face.size()]);
    }
    if(find(adj_edge_of_vertex[new_face[i]].begin(), adj_edge_of_vertex[new_face[i]].end(), new_edge_of_face[i]) ==
        adj_edge_of_vertex[new_face[i]].end()) {
      adj_edge_of_vertex[new_face[i]].push_back(new_edge_of_face[i]);
    }
  }
  for(size_t i = 0; i < new_edge_of_face.size(); ++i) {
    if(find(adj_face_of_edge[new_edge_of_face[i]].begin(), adj_face_of_edge[new_edge_of_face[i]].end(), new_face_index) ==
        adj_face_of_edge[new_edge_of_face[i]].end()) {
      adj_face_of_edge[new_edge_of_face[i]].push_back(new_face_index);
    }
  }
}

int dzw::model::mesh_model::add_face(int edge_index, int vertex_index) {
  std::vector<int> new_face(3);
  new_face[0] = vertex_of_edge[edge_index][0];
  new_face[1] = vertex_of_edge[edge_index][1];
  new_face[2] = vertex_index;
  ///add a new face.
  vertex_of_face.push_back(new_face);
  face_status.push_back(true);
  int new_face_index = vertex_of_face.size() - 1;
  ///add new edges.
  std::vector<int> new_edge_of_face(3);
  for(size_t i = 0; i < new_face.size(); ++i) {
    int new_edge_index = find_edge_index(new_face[i], new_face[(i + 1) % new_face.size()]);
    new_edge_of_face[i] = new_edge_index;
  }
  ///add edge index for new face.
  edge_of_face.push_back(new_edge_of_face);
  ///compute face normal and face area.
  std::vector<matrix<double> > triangle_coord(vertex_of_face[new_face_index].size());
  for(size_t j = 0; j < triangle_coord.size(); ++j) {
    triangle_coord[j] = vertex_coord[vertex_of_face[new_face_index][j]];
  }
  zjucad::matrix::matrix<double> new_face_normal;
  dzw::common::normal(triangle_coord, new_face_normal);
  face_normal.push_back(new_face_normal);
  double new_face_area = dzw::common::area(triangle_coord);
  face_area.push_back(new_face_area);
  ///update adjacent information.
  for(size_t i = 0; i < new_face.size(); ++i) {
    adj_face_of_vertex[new_face[i]].push_back(new_face_index);
  }
  for(size_t i = 0; i < new_face.size(); ++i) {
    if(find(adj_edge_of_vertex[new_face[i]].begin(), adj_edge_of_vertex[new_face[i]].end(), new_edge_of_face[(i + 2) % new_face.size()]) ==
        adj_edge_of_vertex[new_face[i]].end()) {
      adj_edge_of_vertex[new_face[i]].push_back(new_edge_of_face[(i + 2) % new_face.size()]);
    }
    if(find(adj_edge_of_vertex[new_face[i]].begin(), adj_edge_of_vertex[new_face[i]].end(), new_edge_of_face[i]) ==
        adj_edge_of_vertex[new_face[i]].end()) {
      adj_edge_of_vertex[new_face[i]].push_back(new_edge_of_face[i]);
    }
  }
  for(size_t i = 0; i < new_edge_of_face.size(); ++i) {
    if(find(adj_face_of_edge[new_edge_of_face[i]].begin(), adj_face_of_edge[new_edge_of_face[i]].end(), new_face_index) ==
        adj_face_of_edge[new_edge_of_face[i]].end()) {
      adj_face_of_edge[new_edge_of_face[i]].push_back(new_face_index);
    }
  }
}

int dzw::model::mesh_model::delete_vertex(int vertex_index) {
  vertex_status[vertex_index] = false;
  cout << vertex_index << " : " << adj_face_of_vertex[vertex_index].size() << std::endl;
  std::vector<int> delete_face_list, delete_edge_list;
  for(size_t i = 0; i < adj_face_of_vertex[vertex_index].size(); ++i) {
    delete_face_list.push_back(adj_face_of_vertex[vertex_index][i]);
  }
  for(size_t i = 0; i < adj_edge_of_vertex[vertex_index].size(); ++i) {
    delete_edge_list.push_back(adj_edge_of_vertex[vertex_index][i]);
  }
  for(size_t i = 0; i < delete_face_list.size(); ++i) {
    delete_face(delete_face_list[i]);
  }
  for(size_t i = 0; i < delete_edge_list.size(); ++i) {
    delete_edge(delete_edge_list[i]);
  }
  return 0;
}

int dzw::model::mesh_model::delete_face(int face_index) {
  face_status[face_index] = false;
  for(size_t i = 0; i < vertex_of_face[face_index].size(); ++i) {
    int vertex_index = vertex_of_face[face_index][i];
    adj_face_of_vertex[vertex_index].erase(find(adj_face_of_vertex[vertex_index].begin(), adj_face_of_vertex[vertex_index].end(), face_index));
  }
  for(size_t i = 0; i < edge_of_face[face_index].size(); ++i) {
    int edge_index = edge_of_face[face_index][i];
    adj_face_of_edge[edge_index].erase(find(adj_face_of_edge[edge_index].begin(), adj_face_of_edge[edge_index].end(), face_index));
  }
  return 0;
}

int dzw::model::mesh_model::delete_edge(int edge_index) {
  edge_status[edge_index] = false;
  for(size_t i = 0; i < vertex_of_edge[edge_index].size(); ++i) {
    int vertex_index = vertex_of_edge[edge_index][i];
    adj_edge_of_vertex[vertex_index].erase(find(adj_edge_of_vertex[vertex_index].begin(), adj_edge_of_vertex[vertex_index].end(), edge_index));
  }
  return 0;
}

int dzw::model::mesh_model::update_vertex_coord(int index, const zjucad::matrix::matrix<double>& new_coord) {
  vertex_coord[index] = new_coord;
  return 0;
}

void dzw::model::mesh_model::set_one_feature_line(const std::vector<int>& one_feature_line) {
  assert(one_feature_line.size() > 1);
  for(size_t i = 1; i < one_feature_line.size(); ++i) {
    int v1 = one_feature_line[i - 1], v2 = one_feature_line[i];
    assert(v1 < vertex_coord.size() && v2 < vertex_coord.size());
    int edge_index = get_edge_index(v1, v2);
    assert(edge_index != -1);
    feature_vertex[v1] = feature_vertex[v2] = true;
    feature_edge[edge_index] = true;
  }
  feature_line_vec.push_back(one_feature_line);
}

void dzw::model::mesh_model::set_feature_edge(int v1, int v2) {
  assert(v1 < vertex_coord.size() && v2 < vertex_coord.size());
  int edge_index = get_edge_index(v1, v2);
  assert(edge_index != -1);
  feature_vertex[v1] = feature_vertex[v2] = true;
  feature_edge[edge_index] = true;
}
