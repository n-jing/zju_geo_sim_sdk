#include "OpenFlipper.h"

#include <float.h> // for FLT_MIN

#include <memory>
#include <vector>

#include "OpenMesh/src/OpenMesh/Core/IO/MeshIO.hh"

#include "OpenFlipper/common/bsp/TriangleBSPT.hh" // for IsotropicRemesherT.hh
#include "ObjectTypes/TriangleMesh/TriangleMeshTypes.hh" // for TriMesh
#include "Plugin-IsotropicRemesher/IsotropicRemesherT.hh"
#include "OpenMesh/Core/IO/importer/ImporterT.hh"

using namespace std;
using namespace OpenMesh;

int convert(const tri_mesh &src, TriMesh &dst)
{
	IO::ImporterT<TriMesh> imp(dst);
	imp.prepare();
	dst.clear();
	vector<TriMesh::VertexHandle> vh(src.vert_num);
	for(int vi = 0; vi < src.vert_num; ++vi) {
		vh[vi] = dst.add_vertex(TriMesh::Point(src.vertex[3*vi+0],
										   src.vertex[3*vi+1],
										   src.vertex[3*vi+2]));
	}
	for(int fi = 0; fi < src.face_num; ++fi) {
		dst.add_face(vh[src.face[3*fi+0]],
					 vh[src.face[3*fi+1]],
					 vh[src.face[3*fi+2]]);
	}
	imp.finish();

    dst.request_vertex_normals();
    dst.request_face_normals();
    dst.request_vertex_status();
    dst.request_face_status();
    dst.request_edge_status();
	return 0;
}

int convert(const TriMesh &src, tri_mesh &dst)
{
	for(unsigned int vi = 0; vi < src.n_vertices(); ++vi) {
		TriMesh::Point v = src.point(VertexHandle(vi));
		for(int d = 0; d < 3; ++d)
			dst.vertex[vi*3+d] = v[d];
	}
	for(unsigned int fi = 0; fi < src.n_faces(); ++fi) {
		TriMesh::CFVIter fv_it=src.cfv_iter(FaceHandle(fi));
		for (int j = 0; fv_it; ++fv_it, ++j)
			dst.face[fi*3+j] = fv_it.handle().idx();
	}
	return 0;
}

int IsotropicRemesher0(
	const tri_mesh *input,
	double edge_length,
	int *vert_num, int *face_num,
	void **ctx)
{
	IsotropicRemesher< TriMesh > remesher;
	auto_ptr<TriMesh> mesh(new TriMesh);
	if(convert(*input, *mesh))
		return __LINE__;

	remesher.remesh(*mesh, edge_length);

	*vert_num = mesh->n_vertices();
	*face_num = mesh->n_faces();
	*ctx = mesh.release();
	return 0;
}

int IsotropicRemesher1(
	void *ctx,
	tri_mesh *output)
{
	auto_ptr<TriMesh> mesh(reinterpret_cast<TriMesh *>(ctx));
	if(convert(*mesh, *output))
		return __LINE__;
	return 0;
}
