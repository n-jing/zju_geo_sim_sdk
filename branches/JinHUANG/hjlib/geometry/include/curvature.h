#ifndef HJ_GEOMETRY_CURVATURE_H_
#define HJ_GEOMETRY_CURVATURE_H_

#include "conf.h"

namespace hj { namespace geometry {

//! @brief standard normal cycle method [Alliez et al., Anisotropic
//! Polygonal Remeshing, 2003].

//! @param nuv normal, max and min principle directions in right handl
//! orientation.

//! @param cuv normal, max and min principle value from the eigen
//! problem.

//! @note This method under-estimates the magnitude, may be not as
//! good as fitting strategy.

int curvature_by_normal_cycle(const matrixr_t &pts, const matrixs_t &trs,
                              matrixr_t &nuv, matrixr_t &curv);

}}

#endif
