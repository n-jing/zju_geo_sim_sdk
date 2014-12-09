/**   
* @file feature_line.cpp
* @brief TODO
* @author dangzw
* @date Oct 22, 2011 3:32:58 PM 
* @version V1.0   
*/

#include <fstream>
#include <limits>
#include <queue>
#include <limits>
#include <sstream>

#include "../Common/MyException.h"
#include "FeatureLine.h"
#include "MeshModel.h"

using namespace std;
using namespace dzw::common;
using namespace zjucad::matrix;


dzw::model::feature_line::feature_line(const boost::property_tree::ptree& opt_) : opts(opt_){
  if(zjucad::has("input/fl_file_name.value", opts)){
    std::string file_name = opts.get<string>("input/fl_file_name.value");
    load_feature_line(file_name);
  }
}

int dzw::model::feature_line::load_feature_line(const std::string & file_name){
  ifstream ifs(file_name.c_str());
  if(!ifs)
    throw my_exception("Can not open the fl model.");
  int nVertex = 0, nFace = 0; /// Vertex number and face number.
  char format[10];
  std::string str;
  while(getline(ifs, str)) {
    if(str[0] == '#')   continue;
    if(str[0] == 'g')   continue;
    if(str[0] == 'v' && str[1] == ' ') {
      ++nVertex;
      continue;
    }
  }
  ifs.clear();
  ifs.seekg(0, std::ifstream::beg);

  /// Store vertex information
  feature_point_coord.resize(nVertex);
  /// Store edge information
  size_t vn = 0;
  while(getline(ifs, str)) {
    if(str[0] == '#')   continue;
    if(str[0] == 'g')   continue;
    if(str[0] == 'o')   continue;
    // read the vertex information here.
    if(str[0] == 'v' && str[1] == ' ') {
      feature_point_coord[vn].resize(3, 1);
      sscanf(str.c_str(), "%s %lf %lf %lf", format, &feature_point_coord[vn][0], &feature_point_coord[vn][1], &feature_point_coord[vn][2]);
      ++vn;
      continue;
    }
    // read the face information here.
    if((str[0] == 'f' || str[0] == 'l') && str[1] == ' ') {
      std::stringstream ss;
      ss << str;
      string temp;
      ss >> temp;
      std::cout << temp << std::endl;
      int index_1 = 0, index_2;
      ss >> index_1;
      while(ss >> index_2) {
        std::vector<int> edge_i(2);
        edge_i[0] = index_1 - 1; edge_i[1] = index_2 - 1;
        feature_edge_vec.push_back(edge_i);
        index_1 = index_2;
      }
    }
  }
  ifs.close();
  std::cout << "feature edge size : " << feature_edge_vec.size() << std::endl;
  compute_bounding_box();
  std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  find_adjacent_info();
  std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  return 0;
}

void dzw::model::feature_line::compute_bounding_box() {
  model_bounding_center.resize(3, 1);
  matrix<double> max_point = ones<double>(3, 1) * std::numeric_limits<double>::infinity() * -1.0;
  matrix<double> min_point = ones<double>(3, 1) * std::numeric_limits<double>::infinity();
  for(size_t i = 0; i < 3; ++i){
    for(size_t j = 0; j < feature_point_coord.size(); ++j) {
      if(max_point(i, 0) < feature_point_coord[j][i])
        max_point(i, 0) = feature_point_coord[j][i];
      if(min_point(i, 0) > feature_point_coord[j][i])
        min_point(i, 0) = feature_point_coord[j][i];
    }
  }
  model_bounding_center = (max_point + min_point) / 2;
  model_bounding_radius = norm(max_point - min_point) / 2;
}

void dzw::model::feature_line::find_adjacent_info() {
  adj_edge_index_of_vertex.resize(feature_point_coord.size());
  for(size_t i = 0; i < feature_edge_vec.size(); ++i) {
    adj_edge_index_of_vertex[feature_edge_vec[i][0]].push_back(i);
    adj_edge_index_of_vertex[feature_edge_vec[i][1]].push_back(i);
  }
}

int dzw::model::feature_line::get_edge_index(int vertex_1, int vertex_2) const {
  const std::vector<int>& adj_edge_index_1 = adj_edge_index_of_vertex[vertex_1];
  const std::vector<int>& adj_edge_index_2 = adj_edge_index_of_vertex[vertex_2];
  for(size_t i = 0; i < adj_edge_index_1.size(); ++i) {
    for(size_t j = 0; j < adj_edge_index_2.size(); ++j) {
      if(adj_edge_index_1[i] == adj_edge_index_2[j])
        return adj_edge_index_1[i];
    }
  }
  return -1;
}
