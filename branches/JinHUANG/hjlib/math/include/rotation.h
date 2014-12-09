#ifndef HJ_ROTATION_H_
#define HJ_ROTATION_H_

#include "conf.h"

#include <zjucad/matrix/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

// axis is a normalized length 3 vector
// rot is column major 3x3 matrix
// output should be pre-allocated
void HJ_MATH_API axis_angle_to_rot(const double *axis, double angle, double *rot);

/**
   @param axisangle will be normalized

   @return 1 means too small rotation (angle*angle < eps)
 */
int HJ_MATH_API axisangle_to_rot(double *axisangle, double *rot, double eps);

/**
   @param rot will not be modified
 */
void HJ_MATH_API rot_to_axis_angle_Rodrigiues(const double *rot, double *axis, double *angle, double eps);

void HJ_MATH_API rot_to_axis_angle_direct(const double *rot, double *axis, double *angle, double eps);

void HJ_MATH_API M33_to_Euler_XYZ(const double *R, double *XYZ);
void HJ_MATH_API Euler_XYZ_to_M33(const double *XYZ, double *R);

#ifdef __cplusplus
}
#endif

#endif
