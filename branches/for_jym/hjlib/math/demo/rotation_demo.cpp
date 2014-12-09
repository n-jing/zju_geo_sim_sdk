#include <iostream>
using namespace std;


#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
using namespace zjucad::matrix;

#include "rotation.h"

int main(int argc, char *argv[])
{
	matrix<double> axis = rand<double>(3, 1);
	axis /= norm(axis);
	double angle = 1.2;
	cout << axis << angle << endl;

	matrix<double> R(3, 3);

	axis_angle_to_rot(&axis[0], angle, &R[0]);
	cout << R << norm(R*trans(R)-eye<double>(3)) << endl;;

	rot_to_axis_angle_Rodrigiues(&R[0], &axis[0], &angle, 1e-10);
	cout << axis << angle << endl;

	rot_to_axis_angle_direct(&R[0], &axis[0], &angle, 1e-10);
	cout << axis << angle << endl;

	double Ra[] = {
		1, 3.33067e-016, 6.66134e-016,
		0, 1, 2.66454e-015,
		-4.44089e-016, -1.77636e-015, 1};
	rot_to_axis_angle_direct(Ra, &axis[0], &angle, 1e-20);
	cout << "check error " << axis << angle << endl;
	rot_to_axis_angle_Rodrigiues(Ra, &axis[0], &angle, 1e-10);
	cout << axis << angle << endl;

	return 0;
}
