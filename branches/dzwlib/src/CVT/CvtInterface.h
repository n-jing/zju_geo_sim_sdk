/**   
* @file CvtInterface.h 
* @brief TODO
* @author dangzw
* @date May 25, 2012 5:22:15 PM 
* @version V1.0   
*/

#ifndef CVTINTERFACE_H_
#define CVTINTERFACE_H_

#include "../Common/BaseData.h"

/** Not use
  */
int get_rvd_mesh(const std::string& file_name,
    const dzw::common::VertexCoordMat& seeds_coord,
    dzw::common::VertexCoordMat& rvd_vertex,
    dzw::common::VertexOfFaceMat& rvd_face,
    std::vector<std::vector<int> >& rvd_classify,
    dzw::common::VertexCoordMat& rdt_vertex,
    dzw::common::VertexOfFaceMat& rdt_face);
/** Not use
  */
int save_rdt_mesh(const std::string& base_file_name,
    const std::string& rvd_file_name,
    const std::string& rdt_file_name,
    const dzw::common::VertexCoordMat& seeds_coord);

/** Get the value and gradient of CVT energy function.
  * @param file_name the file name of original mesh
  * @param seeds_coord the seed vertices coordinate
  * @param gradient the gradient of CVT energy function
  * @return the value of CVT energy function
  */

double get_rvd_f_g(const std::string& file_name,
    const dzw::common::VertexCoordMat& seeds_coord,
    std::vector<double>& gradient);

//double get_rvd_f_g(const std::string& file_name,
//    const std::string& pts_file);

#endif /* CVTINTERFACE_H_ */
