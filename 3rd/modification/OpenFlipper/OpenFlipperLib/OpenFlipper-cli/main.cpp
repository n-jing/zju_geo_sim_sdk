#include <iostream>
#include <vector>
#include <iterator>

#include "../OpenFlipperLib/OpenFlipper.h"

using namespace std;

int main(int argc, char *argv[])
{
	double vertex[9] = {
		0, 0, 0, 1, 0, 0, 0, 1, 0
	};
	int face[3] = {
		0, 1, 2
	};
	tri_mesh input;
	input.vertex = vertex;
	input.face = face;
	input.vert_num = 3;
	input.face_num = 1;

	int vert_num, face_num;
	void *ctx;
	if(IsotropicRemesher0(&input, 0.2, &vert_num, &face_num, &ctx))
		return __LINE__;
	vector<double> vertex1(vert_num*3);
	vector<int> face1(face_num*3);
	tri_mesh output;
	output.vertex = &vertex1[0];
	output.face = &face1[0];
	IsotropicRemesher1(ctx, &output);

	copy(vertex1.begin(), vertex1.end(), ostream_iterator<double>(cout, " "));
	cout << endl;
	copy(face1.begin(), face1.end(), ostream_iterator<int>(cout, " "));
	cout << endl;
	cerr << "success." << endl;
	
	return 0;
}
