/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *      Copyright (C) 2001-2010 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openflipper.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenFlipper.                                        *
 *                                                                           *
 *  OpenFlipper is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenFlipper is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenFlipper. If not,                                  *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 9595 $                                                       *
 *   $Author: moebius $                                                      *
 *   $Date: 2010-06-17 12:48:23 +0200 (Do, 17. Jun 2010) $                   *
 *                                                                           *
\*===========================================================================*/



//=============================================================================
//
//  CLASS GLMatrixT
//
//=============================================================================


#ifndef ACG_GLMATRIX_HH
#define ACG_GLMATRIX_HH


//== INCLUDES =================================================================


#include "Matrix4x4T.hh"
#include "../Config/ACGDefines.hh"
#include <math.h>


namespace ACG {

	      
//== CLASS DEFINITION =========================================================


/** enum to chose whether to multiply new matrices from right
    (default) or left to the current matrix */
enum MultiplyFrom { MULT_FROM_RIGHT, MULT_FROM_LEFT };



/// 4x4 matrix implementing OpenGL commands.
template <class Scalar>
class GLMatrixT : public Matrix4x4T<Scalar>
{
public:
   
  typedef VectorT<Scalar, 3> Vec3;


  /// constructor: uninitialized values
  GLMatrixT() {}
 
  /// construct from other matrix type
  template <class OtherScalar>
  inline GLMatrixT(const GLMatrixT<OtherScalar> _rhs) 
    : Matrix4x4T<Scalar>(_rhs)
  {}

  /** setup matrix using an array of N*N scalar values.
      elements are ordered 'column first' (like OpenGL) */
  inline GLMatrixT(const Scalar _array[16]) : Matrix4x4T<Scalar>(_array) {}

  /// destructor
  ~GLMatrixT() {}


  /// assignement from other matrix type
  template<typename otherScalar>
  inline GLMatrixT<Scalar>& operator=(const GLMatrixT<otherScalar>& _rhs)
  { Matrix4x4T<Scalar>::operator=(_rhs); return *this; }

  /// assignement from other matrix type
  template<typename otherScalar>
  inline GLMatrixT<Scalar>& operator=(const Matrix4x4T<otherScalar>& _rhs)
  { Matrix4x4T<Scalar>::operator=(_rhs); return *this; }


      
  /// multiply self with scaling matrix (x,y,z)
  inline void scale( Scalar _x, Scalar _y, Scalar _z, 
		     MultiplyFrom _mult_from = MULT_FROM_RIGHT );
  /// multiply self with scaling matrix (x,y,z)
  inline void scale( const Vec3& _v,
		     MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    scale(_v[0], _v[1], _v[2], _mult_from);
  }


  /// multiply self with translation matrix (x,y,z)
  inline void translate( Scalar _x, Scalar _y, Scalar _z,
			 MultiplyFrom _mult_from = MULT_FROM_RIGHT );
  /// multiply self with translation matrix (x,y,z)
  inline void translate( const Vec3& _v,
			 MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    translate(_v[0], _v[1], _v[2], _mult_from);
  }


  /** multiply self with a rotation matrix
      (angle in degree, arbitrary axis given by xyz) */
  void rotate( Scalar angle, Scalar x, Scalar y, Scalar z,
	       MultiplyFrom _mult_from = MULT_FROM_RIGHT );
  /** multiply self with a rotation matrix
      (angle in degree, arbitrary axis given by xyz) */
  void rotate( Scalar _angle, const Vec3& _axis,
	       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotate(_angle, _axis[0], _axis[1], _axis[2], _mult_from);
  }



  /// multiply self with a rotation matrix (angle in degree, x-axis) 
  inline void rotateX( Scalar _angle,
		       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotateXYZ( X, _angle, _mult_from );
  }

  /// multiply self with a rotation matrix (angle in degree, y-axis) 
  inline void rotateY( Scalar _angle,
		       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotateXYZ( Y, _angle, _mult_from );
  }

  /// multiply self with a rotation matrix (angle in degree, z-axis) 
  inline void rotateZ( Scalar _angle,
		       MultiplyFrom _mult_from = MULT_FROM_RIGHT ) {
    rotateXYZ( Z, _angle, _mult_from );
  }




  /** multiply self with a viewing transformation given by 
      eye position, reference point (center) and an up vector.
      (similar to gluLookAt) */
  void lookAt(const Vec3& eye,
	      const Vec3& center,
	      const Vec3& up);

  /// multiply self from left with inverse lookAt matrix
  void inverse_lookAt(const Vec3& eye,
		      const Vec3& center,
		      const Vec3& up);


  /// multiply self with a perspective projection matrix 
  void perspective(Scalar fovY, Scalar aspect, 
		   Scalar near_plane, Scalar far_plane);

  /// multiply self from left with inverse of perspective projection matrix 
  void inverse_perspective(Scalar fovY, Scalar aspect,
			   Scalar near_plane,Scalar far_plane);

  /// multiply self with a perspective projection matrix 
  void frustum(Scalar left, Scalar right,
	       Scalar bottom, Scalar top,
	       Scalar near_plane, Scalar far_plane);

  /// multiply self from left with inverse of perspective projection matrix 
  void inverse_frustum(Scalar left,Scalar right,
		       Scalar bottom, Scalar top,
		       Scalar near_plane, Scalar far_plane);

  /// multiply self with orthographic projection matrix
  void ortho(Scalar left, Scalar right,
	     Scalar bottom, Scalar top,
	     Scalar near_plane, Scalar far_plane);

  /// multiply self from left with inverse orthographic projection matrix
  void inverse_ortho(Scalar left, Scalar right,
		     Scalar bottom, Scalar top,
		     Scalar near_plane, Scalar far_plane);




  //----------------------------------------------------- overloaded operators 

  GLMatrixT& operator+= ( const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::operator+=(_rhs); return *this;
  }
  GLMatrixT& operator-= ( const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::operator-=(_rhs); return *this;
  }
  GLMatrixT& operator*= ( const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::operator*=(_rhs); return *this;
  }
  GLMatrixT& leftMult(const Matrix4x4T<Scalar>& _rhs) {
    Matrix4x4T<Scalar>::leftMult(_rhs); return *this;
  }

  GLMatrixT operator+ (const Matrix4x4T<Scalar>& _rhs) const {
    return GLMatrixT<Scalar>(*this) += _rhs;
  }    
  GLMatrixT operator- (const Matrix4x4T<Scalar>& _rhs) const {
    return GLMatrixT<Scalar>(*this) -= _rhs;
  }
  GLMatrixT operator*(const Matrix4x4T<Scalar>& _rhs) const {
    return GLMatrixT<Scalar>(*this) *= _rhs;
  }

  template <typename T>
  inline VectorT<T,4> operator*(const VectorT<T,4>& _v) const {
    return Matrix4x4T<Scalar>::operator*(_v);
  }



private:

  enum Axis { X, Y, Z };
  void rotateXYZ( Axis _axis, Scalar _angle, MultiplyFrom _mult_from );
};




//=============================================================================


/// typedef
typedef GLMatrixT<float>  GLMatrixf;
/// typedef
typedef GLMatrixT<double> GLMatrixd;


//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_GLMATRIX_C)
#define ACG_GLMATRIX_TEMPLATES
#include "GLMatrixT.cc"
#endif
//=============================================================================
#endif // ACG_GLMATRIX_HH defined
//=============================================================================

