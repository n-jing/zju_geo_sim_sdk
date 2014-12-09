/**   
* @file BaseData.cpp 
* @brief TODO
* @author dangzw
* @date Nov 19, 2011 7:58:30 PM 
* @version V1.0   
*/

#include "BaseData.h"

//int N_sym = 4;

using namespace zjucad::matrix;

double dzw::common::mod(double v1, double v2){
  assert(fabs(v2) > nearly_zero && v2 > 0.0);
  int c = v1 / v2;
  double temp = v1 - v2 * double(c);
  if(temp < -0.5 * v2)
    return temp + 0.5 * v2;
  else if(temp > 0.5 * v2){
    return temp - 0.5 * v2;
  }
  else
    return temp;
}

double dzw::common::angle(const zjucad::matrix::matrix<double>& v1, const zjucad::matrix::matrix<double>& v2) {
  double cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
  cos_angle = cos_angle > 1.0 ? 1.0 : cos_angle;
  cos_angle = cos_angle < -1.0 ? -1.0 : cos_angle;
  return acos(cos_angle);
}

double dzw::common::distance_vertex_to_edge(const zjucad::matrix::matrix<double>& v1,
    const zjucad::matrix::matrix<double>& v2,
    const zjucad::matrix::matrix<double>& v,
    zjucad::matrix::matrix<double>& intersect_v) {
  /*
  double A1, A2, B1, B2, C1, C2, D1, D2;
  double x1 = v1(0, 0), y1 = v1(1, 0), z1 = v1(2, 0);
  double x2 = v2(0, 0), y2 = v2(1, 0), z2 = v2(2, 0);
  A1 = y2 - y1; B1 = x1 - x2; C1 = 0.0;
  D1 = x2 * y1 - x1 * y2;
  A2 = z2 - z1; B2 = 0.0; C2 = x1 - x2;
  D2 = x2 * z1 - x1 * z2;
  zjucad::matrix::matrix<double> n1(3, 1), n2(3, 1);
  n1(0, 0) = A1; n1(1, 0) = B1; n1(2, 0) = C1;
  n2(0, 0) = A2; n2(1, 0) = B2; n2(2, 0) = C2;
  double d1 = A1 * v(0, 0) + B1 * v(1, 0) + C1 * v(2, 0) + D1;
  double d2 = A2 * v(0, 0) + B2 * v(1, 0) + C2 * v(2, 0) + D2;
  double dis = norm(d1 * n2 - d2 * n1) * norm(cross(n1, n2));
  return dis;*/
  double param_t = dot(v1 - v, v1 - v2) / dot(v1 - v2, v1 - v2);
  if(param_t >= 0.0 && param_t <= 1.0)
    intersect_v = v1 + param_t * (v2 - v1);
  else if(param_t < 0.0)
    intersect_v = v1;
  else
    intersect_v = v2;
  return norm(v - intersect_v);
}

double dzw::common::area(const std::vector<zjucad::matrix::matrix<double> >& triangle) {
  /// area = |edge1 X edge2| / 2;
  assert(triangle.size() == 3);
  return 0.5 * norm(cross(triangle[1] - triangle[0],
      triangle[2] - triangle[1]));
}

int dzw::common::normal(const std::vector<zjucad::matrix::matrix<double> >& triangle, zjucad::matrix::matrix<double>& normal_) {
  ///Cross two edge vector of triangle to compute the face normal.
  assert(triangle.size() == 3);
  normal_ = cross(triangle[1] - triangle[0],
      triangle[2] - triangle[1]);
  if(norm(normal_) < nearly_zero){
    normal_ = zeros<double>(3, 1);
    return 1;
  }
  ///Normalize the face normal.
  normal_ /= norm(normal_);
  return 0;
}

int dzw::common::get_plane_equation(const std::vector<zjucad::matrix::matrix<double> >& vertex_vec, zjucad::matrix::matrix<double>& plane_nomal, double& plane_d) {
  plane_nomal = cross(temp(vertex_vec[2] - vertex_vec[1]),temp(vertex_vec[0] - vertex_vec[1]));
  plane_nomal /= norm(plane_nomal);
  plane_d = 0.0;
  for(size_t i = 0;i < plane_nomal.size(1); ++i){
    plane_d -= (plane_nomal(i,0) * vertex_vec[0](i,0));
  }
  return 0;
}

int dzw::common::get_pedal_vertex_to_plane(const zjucad::matrix::matrix<double>& out_vertex,
    const zjucad::matrix::matrix<double>& plane_normal,
    double plane_d, zjucad::matrix::matrix<double>& pedal) {
  double param_t = (-plane_d - dot(plane_normal, out_vertex)) / dot(plane_normal, plane_normal);
  pedal = out_vertex + param_t * plane_normal;
  return 0;
}

bool dzw::common::is_vertex_in_triangle(const zjucad::matrix::matrix<double>& vertex_coord,
    const std::vector<zjucad::matrix::matrix<double> >& triangle_coord) {
  std::vector<zjucad::matrix::matrix<double> > norm_vec(triangle_coord.size());
  for(size_t i = 0;i < triangle_coord.size(); ++i){
    norm_vec[i] = triangle_coord[i] - vertex_coord;
  }
  zjucad::matrix::matrix<double> v0 = cross(norm_vec[0],norm_vec[1]);
  if(norm(v0) != 0){
    for(size_t i = 1;i < norm_vec.size(); ++i){
      zjucad::matrix::matrix<double> vi = cross(norm_vec[i],norm_vec[(i + 1) % norm_vec.size()]);
      if(dot(v0,vi) < 0)
        return false;
    }
  }
  else{
    if(dot(norm_vec[0],norm_vec[1]) > 0)
      return false;
  }
  return true;
}

bool dzw::common::is_vertex_in_polygon(const zjucad::matrix::matrix<double>& vertex_coord,
                          const zjucad::matrix::matrix<double>& polygon_coord)
{
  for(size_t i = 1; i < (polygon_coord.size(2) - 1); ++i) {
    std::vector<zjucad::matrix::matrix<double> > triangle_coord(3);
    triangle_coord[0] = polygon_coord(colon(), 0);
    triangle_coord[1] = polygon_coord(colon(), i);
    triangle_coord[2] = polygon_coord(colon(), i + 1);
    if(is_vertex_in_triangle(vertex_coord, triangle_coord))
      return true;
  }
  return false;
}

bool dzw::common::is_intersect_line_polygon(const zjucad::matrix::matrix<double>& origin,
                           const zjucad::matrix::matrix<double>& des,
                           const zjucad::matrix::matrix<double>& polygon_coord)
{
  zjucad::matrix::matrix<double> direction = des - origin;
  direction /= norm(direction);
  zjucad::matrix::matrix<double> e1 = polygon_coord(colon(), 0) - polygon_coord(colon(), 1);
  zjucad::matrix::matrix<double> e2 = polygon_coord(colon(), 2) - polygon_coord(colon(), 1);
  zjucad::matrix::matrix<double> p_norm = cross(e1, e2);
  double parrell = dot(p_norm,direction);
  if(fabs(parrell) != dzw::common::nearly_zero){
    double ray_t = dot(p_norm, polygon_coord(colon(), 1) - origin) / parrell;
    if(ray_t < 0 || ray_t > norm(des - origin))
      return false;
    zjucad::matrix::matrix<double> itersect_point = ray_t * direction + origin;
    if(is_vertex_in_polygon(itersect_point, polygon_coord))
      return true;
  }
  return false;
}

bool dzw::common::is_intersect_two_polygons(const zjucad::matrix::matrix<double>& polygon_coord_1,
                           const zjucad::matrix::matrix<double>& polygon_coord_2)
{
  for(size_t i = 1; i < polygon_coord_1.size(2); ++i) {
    if(is_intersect_line_polygon(polygon_coord_1(colon(), i - 1), polygon_coord_1(colon(), i), polygon_coord_2))
      return true;
  }
  for(size_t i = 1; i < polygon_coord_2.size(2); ++i) {
    if(is_intersect_line_polygon(polygon_coord_2(colon(), i - 1), polygon_coord_2(colon(), i), polygon_coord_1))
      return true;
  }
  return false;
}
