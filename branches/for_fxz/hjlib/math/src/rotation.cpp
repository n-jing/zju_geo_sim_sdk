#include "rotation.h"

#include <algorithm>
using namespace std;

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
using namespace zjucad::matrix;

template <typename T>
static inline void cross_matrix(const T &u,
				  matrix<double> &m)
{
	m[0] = m[4] = m[8] = 0;
	m[1] = u[2];
	m[2] = -u[1];
	m[3] = -u[2];
	m[5] = u[0];
	m[6] = u[1];
	m[7] = -u[0];
}

void axis_angle_to_rot(const double *axis, double angle, double *rot)
{
  matrix<double> A(3, 3);
  itr_matrix<double *> R(3, 3, rot);
	const double c = cos(angle), s = sin(angle);
	cross_matrix(axis, A);
	R = eye<double>(3)+s*A+(1-c)*(A*A);
}

int axisangle_to_rot(double *axisangle, double *rot, double eps)
{
	double angle = 0;
	for(int i = 0; i < 3; ++i)
		angle += axisangle[i]*axisangle[i];
	if(angle < eps) {
		fill(rot, rot+9, 0);
		rot[0] = rot[4] = rot[8] = 1;
		return 1;
	}
	else {
		angle = sqrt(angle);
		for(int i = 0; i < 3; ++i)
			axisangle[i] /= angle;
	}
	axis_angle_to_rot(axisangle, angle, rot);
	return 0;
}

// Has problem when near pi
void rot_to_axis_angle_Rodrigiues(const double *rot, double *axis, double *angle, double eps)
{
	const itr_matrix<const double *> R(3, 3, rot);
	matrix<double> A = (R-trans(R))/2;
	itr_matrix<double *> w(3, 1, axis);
	w[0] = A[5];
	w[1] = A[6];
	w[2] = A[1];

	*angle = norm(w);
	if(*angle < eps) {
		*angle = 0;
		w[0] = 1;
		w[1] = w[2] = 0;
	}
	else {
		w /= *angle;
		if(*angle > 1)
			*angle = 1;
		*angle = asin(*angle);
	}
}

static inline int sign(double v)
{
	if(v == 0) return 0;
	if(v > 0) return 1;
	if(v < 0) return -1;
}

void rot_to_axis_angle_direct(const double *rot, double *axis, double *angle, double eps)
{
	const itr_matrix<const double *> R(3, 3, rot);
  matrix<double> d(3);
  itr_matrix<double *> w(3, 1, axis);
	for(int i = 0; i < 3; ++i)
		d[i] = R(i, i);
	const double t = sum(d);//1+2cos(angle)

	matrix<double> uR(3), uR2(3);
	uR[0] = R(2,1)-R(1,2);
	uR[1] = R(0,2)-R(2,0);
	uR[2] = R(1,0)-R(0,1);
	for(int i = 0; i < 3; ++i)
		uR2[i] = uR[i]*uR[i];

	const double s2 = sum(uR2);// 4sin^2(angle)
	const double a = (3-t)+s2;

    if(a > eps) {
		*angle = atan2(sqrt(s2), t-1);

		matrix<double> A = R+trans(R)-t*eye<double>(3);
		A = pow(A, 2);
		matrix<double> u2=A*ones<double>(3, 1) + uR2 + 4*d + (1-2*t);

		matrix<double> absu(3);
		matrix<int> signu(3);
		for(int i = 0; i < 3; ++i) {
			absu[i] = sqrt(u2[i]*double(u2[i]>0))/sqrt(a);
			signu[i] = sign(uR[i]);
			if(signu[i] == 0)
				signu[i] = 1; // we don't want the nil case
		}

		for(int i = 0; i < 3; ++i)
			w[i] = absu[i]*signu[i];

		double w2 = dot(w, w);
		if(w2 < eps) {
			w[0] = 1; w[1] = w[2] = 0;
		}
		else
			w=w/sqrt(w2);
	}
	else {
		w[0] = 1; w[1] = w[2] = 0;
        angle=0;
	}
}

void HJ_MATH_API M33_to_Euler_XYZ(const double *R0, double *XYZ)
{
  const itr_matrix<const double *> R(3, 3, R0);
  XYZ[0] = atan2(R(2,1), R(2,2));
  XYZ[1] = atan2(-R(2,0), sqrt(R(2,1)*R(2,1) + R(2,2)*R(2,2)));
  XYZ[2] = atan2(R(1,0), R(0,0));
}

void HJ_MATH_API Euler_XYZ_to_M33(const double *XYZ, double *R0)
{
  itr_matrix<double *> R(3, 3, R0);
  matrix<double> X = eye<double>(3), Y = eye<double>(3), Z = eye<double>(3);

  X(1,1) = cos(XYZ[0]);
  X(1,2) = -sin(XYZ[0]);
  X(2,1) = sin(XYZ[0]);
  X(2,2) = cos(XYZ[0]);

  Y(0,0) = cos(XYZ[1]);
  Y(0,2) = sin(XYZ[1]);
  Y(2,0) = -sin(XYZ[1]);
  Y(2,2) = cos(XYZ[1]);

  Z(0,0) = cos(XYZ[2]);
  Z(0,1) = -sin(XYZ[2]);
  Z(1,0) = sin(XYZ[2]);
  Z(1,1) = cos(XYZ[2]);

  R = Z*temp(Y*X);
}
