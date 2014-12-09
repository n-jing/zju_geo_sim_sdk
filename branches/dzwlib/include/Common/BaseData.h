/**   
* @file BaseData.h
* @brief Define some basic data struct
* @author dangzw
* @date Oct 20, 2011 4:05:20 PM 
* @version V1.0   
*/

#ifndef BASEDATA_H_
#define BASEDATA_H_

/** BaseData of the whole project
 *
 */
#include <zjucad/matrix/matrix.h>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>
#include "VectorOp.h"

namespace dzw{namespace common{
/** @def #define INFO_PRINT 1
 * define a macro, 1 means that program need to print information.
 */
#define INFO_PRINT 1

/**@var const double nearly_zero = 1e-10.
 * If less than nearly_zero, will be regard as zero.
 */
const double nearly_zero = 1e-6;

/**@var cosnt double nearly_zero_angle = 0.035
 * If the angle less than nearly_zero_angle(about 2 degree), it will be regard as zero.
 */
const double nearly_zero_angle = 2 * acos(-1.0) / 180;

/**@var const int N_sys = 4
 * N symmetry rotation field
 */
//extern int N_sym;

/** @enum RETURN_TYPE
 * Enumerate function return information, we can analyze the error type.
 */
enum RETURN_TYPE{
	SUCCESS,          //!< SUCCESS
	OPEN_FILE_ERROR,  //!< OPEN_FILE_ERROR can not open file.
	FILE_FORMAT_ERROR,//!< FILE_FORMAT_ERROR file format is different with what we want.
	NEGATIVE_INDEX,   //!< NEGATIVE_INDEX the index of array is negative.
	DIVIDE_ZERO       //!< DIVIDE_ZERO divided by zero.
};

/** @var typedef zjucad::matrix::matrix<double> PointCoord
 * Store the coordinate of a point.
 * Generally, it is a 3*1 matrix.
 * */
typedef zjucad::matrix::matrix<double> PointCoord;

/** @var typedef zjucad::matrix::matrix<double> VertexCoordMat
 * Store the coordinate of vertex of a model.
 * Generally, it is a 3*n matrix, n is the number of vertex.
 */
typedef std::vector<zjucad::matrix::matrix<double> > VertexCoordMat;

/** @var typedef zjucad::matrix::matrix<int> VertexOfFaceMat
 * Store the vertex index of faces.
 * Generally, it is a 3*n matrix, n is the number of face.
 */
typedef std::vector<std::vector<int> > VertexOfFaceMat;

/** @var typedef std::vector<std::vector<int> > EdgeOfFaceVec
 * Store the edge index of faces.
 * Generally, it is a 3*n matrix, n is the number of face.
 */
typedef std::vector<std::vector<int> > EdgeOfFaceVec;

/** @var typedef std::vector<std::vector<int> > VertexOfEdgeVec
 * Store the vertex index of edges.
 * Generally, it is a 2*n matrix, n is the number of edge.
 */
typedef std::vector<std::vector<int> > VertexOfEdgeVec;

/** @var typedef zjucad::matrix::matrix<double> VertexNormalMat
 * Store the vertex normal.
 * Generally, it is a 3*n matrix, n is the number of vertex.
 */
typedef std::vector<zjucad::matrix::matrix<double> > VertexNormalMat;

/** @var typedef zjucad::matrix::matrix<double> FaceNormalMat
 * Store the face normal.
 * Generally, it is a 3*n matrix, n is the number of face.
 */
typedef std::vector<zjucad::matrix::matrix<double> > FaceNormalMat;

/** @var typedef std::vector<std::vector<int> > AdjFaceOfVertex
 * Store the adjacent face index of vertex.
 * It is likely a linklist.
 */
typedef std::vector<std::vector<int> > AdjFaceOfVertex;


/** @var typedef std::vector<std::vector<int> > AdjFaceOfEdge
 * Store the adjacent face index of edge.
 * It is likely a linklist.
 */
typedef std::vector<std::vector<int> > AdjFaceOfEdge;

/** @var typedef std::vector<std::vector<int> > AdjEdgeOfVertex
 * Store the adjacent edge index of vertex.
 * It is likely a linklist.
 */
typedef std::vector<std::vector<int> > AdjEdgeOfVertex;

/** @var typedef zjucad::matrix::matrix<double> VectorFieldMat
 * Store the surface vector field of model.
 * Generally, it is a 2*n matrix, n is the number of edge, face or vertex.
 */
typedef zjucad::matrix::matrix<double> VectorFieldMat;

/**@struct struct Node
 * Dijkstra node.
 */
struct Node{
  int index;
  double dis;
  Node(int index_, double dis_) : index(index_), dis(dis_){}
};

double mod(double v1, double v2);

//template<typename T>
double angle(const zjucad::matrix::matrix<double>& v1, const zjucad::matrix::matrix<double>& v2);

double distance_vertex_to_edge(const zjucad::matrix::matrix<double>& v1,
    const zjucad::matrix::matrix<double>& v2,
    const zjucad::matrix::matrix<double>& v,
    zjucad::matrix::matrix<double>& intersect_v);

double area(const std::vector<zjucad::matrix::matrix<double> >& triangle);

int normal(const std::vector<zjucad::matrix::matrix<double> >& triangle, zjucad::matrix::matrix<double>& normal_);

int get_plane_equation(const std::vector<zjucad::matrix::matrix<double> >& vertex_vec, zjucad::matrix::matrix<double>& plane_nomal, double& plane_d);

int get_pedal_vertex_to_plane(const zjucad::matrix::matrix<double>& out_vertex,
    const zjucad::matrix::matrix<double>& plane_normal,
    double plane_d, zjucad::matrix::matrix<double>& pedal);

bool is_vertex_in_triangle(const zjucad::matrix::matrix<double>& vertex_coord,
    const std::vector<zjucad::matrix::matrix<double> >& triangle_coord);

bool is_vertex_in_polygon(const zjucad::matrix::matrix<double>& vertex_coord,
                          const zjucad::matrix::matrix<double>& polygon_coord);

bool is_intersect_line_polygon(const zjucad::matrix::matrix<double>& origin,
                           const zjucad::matrix::matrix<double>& des,
                           const zjucad::matrix::matrix<double>& polygon_coord);

bool is_intersect_two_polygons(const zjucad::matrix::matrix<double>& polygon_coord_1,
                           const zjucad::matrix::matrix<double>& polygon_coord_2);

}}

#endif /* BASEDATA_H_ */
