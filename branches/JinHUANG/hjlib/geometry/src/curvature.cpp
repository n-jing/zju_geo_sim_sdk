#include "../include/curvature.h"

#include <iostream>

#include <hjlib/half_edge/half_edge.h>
#include <hjlib/half_edge/container.h>
#include <hjlib/half_edge/builder.h>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

using namespace std;
using namespace zjucad::matrix;

namespace hj { namespace geometry {

class entity_with_id
{
public:
  virtual ~entity_with_id(){}
  operator size_t() const { return id_; }
  operator size_t&() { return id_; }
  size_t id() const { return id_; }
  size_t &id() { return id_; }
protected:
  size_t id_;
};

struct internal_property
{
  struct vert_t : public entity_with_id {
  };
  struct edge_t {
  };
  struct face_t : public entity_with_id {
  };
};

typedef hj::half_edge::half_edge_mesh_t<
  hj::half_edge::std_vector, internal_property> mesh_t;

int curvature_by_normal_cycle(const matrixr_t &pts, const matrixs_t &trs,
                                    matrixr_t &nuv, matrixr_t &curv)
{
  mesh_t m;
  try {
    m = hj::half_edge::build_with_id(trs);
  }
  catch (...) {
    cerr << "# half edge patch fail." << endl;
    return __LINE__;
  }

  const double EPS = 1e-12;

  matrixr_t tri_area(m.faces().size());
  double total_area = 0;
  matrixr_t tri_normal(3, m.faces().size());
  for(size_t fi = 0; fi < trs.size(2); ++fi) {
    const matrixr_t tri = pts(colon(), trs(colon(), fi));
    const matrixr_t e[] = {tri(colon(), 1)-tri(colon(), 0),
                           tri(colon(), 2)-tri(colon(), 0)};
    tri_normal(colon(), fi) = cross(e[0], e[1]);
    const double area2 = norm(tri_normal(colon(), fi));
    if(area2 < EPS) {
      cerr << "# waring degenerated triangle: " << fi << endl;
      continue;
    }
    tri_normal(colon(), fi) /= area2;
    total_area += area2/2;
  }

  // curvature tensor
  matrixr_t ct = zeros<real_t>(3, 3); // unit is 1/len.
  for(auto ei = m.edges().begin(); ei != m.edges().end(); ++ei) {
    if(!ei->face()|| !ei->oppo()->face()) continue; // boundary edge
    const size_t a = ei->face()->id(), b = ei->oppo()->face()->id();
    if(a > b) continue;
    matrixr_t edge = pts(colon(), ei->vert()->id())
      - pts(colon(), ei->oppo()->vert()->id());
    const real_t len = norm(edge);
    if(len < EPS) {
      cerr << "# waring degenerated edge between: " << a << " " << b << endl;
      continue;
    }
    edge /= len;
    double c = dot(tri_normal(colon(), a), tri_normal(colon(), b)),
      s = dot(cross(tri_normal(colon(), a), tri_normal(colon(), b)), edge);
    assert(c*c + s*s - 1 < EPS && c <= 1 && c >= 1);
    if(c > 1) c = 1;
    if(c < -1) c = -1;
    double theta = acos(c);
    if(s < 0)
      theta = -theta;
    const matrixr_t cti = edge*trans(edge);
    ct += cti*(len*theta);
  }
  ct /= total_area;

  nuv = ct;
  curv.resize(3);
  eig(nuv, curv);
  { // order it
    swap(curv[1], curv[2]); //! NOTICE here to order it into nomral, kappa_max, kappa_min
    if(dot(nuv(colon(), 0), tri_normal*ones<real_t>(tri_normal.size(2), 1)) < 0)
      nuv *= -1;
    if(dot(nuv(colon(), 0), cross(nuv(colon(), 1), nuv(colon(), 2))) < 0)
      nuv(colon(), 1) *= -1;
  }
  return 0;
}

}}
