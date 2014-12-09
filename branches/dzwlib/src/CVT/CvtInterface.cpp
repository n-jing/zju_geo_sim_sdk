/**   
* @file CvtInterface.cpp 
* @brief TODO
* @author dangzw
* @date May 28, 2012 11:19:57 AM 
* @version V1.0   
*/

#include <fstream>
#include "common/line_stream.h"
#include "CvtInterface.h"
#include "combinatorics/delaunay.h"
#include "combinatorics/RVD.h"
#include "algebra/F_Lp.h"

using namespace Geex;

class get_primal_triangle {
public:
  get_primal_triangle(std::vector<std::vector<int> >& face_vec_) : face_vec(&face_vec_) {
  }
  void operator()(unsigned int i, unsigned int j, unsigned int k) const {
    std::vector<int> face_i(3);
    face_i[0] = i;
    face_i[1] = j;
    face_i[2] = k;
    (*face_vec).push_back(face_i);
  }

private:
    std::vector<std::vector<int> >* face_vec;
} ;

class get_RVD_facets {
public:
  get_RVD_facets(
      dzw::common::VertexCoordMat& vertex_coord_,
      dzw::common::VertexOfFaceMat& face_vec_,
      std::vector<int>& face_flag_
  ) : vertex_coord(vertex_coord_), face_vec(face_vec_), face_flag(face_flag_),
  cur_v_(0) {
  }
  void operator()(unsigned int iv, Mesh* M) const {
    face_vec.resize(M->nb_facets());
    face_flag.resize(M->nb_facets());
    for(size_t f = 0; f < M->nb_facets(); ++f) {
      face_vec[f].resize(M->facet_end(f) - M->facet_begin(f));
      for(size_t i = M->facet_begin(f); i < M->facet_end(f); ++i) {
        const vec3& v = M->vertex(i) ;
        zjucad::matrix::matrix<double> vertex_i(3, 1);
        for(size_t j = 0; j < vertex_i.size(1); ++j)
          vertex_i[j] = v[j];
        vertex_coord.push_back(vertex_i);
      }
      for(size_t i = M->facet_begin(f); i < M->facet_end(f); ++i) {
        face_vec[f][i] = cur_v_;
        cur_v_++ ;
      }
      face_flag[f] = iv;
    }
  }
private:
  dzw::common::VertexCoordMat& vertex_coord;
  dzw::common::VertexOfFaceMat& face_vec;
  std::vector<int>& face_flag;
  mutable unsigned int cur_v_ ;
} ;

void set_pts(const dzw::common::VertexCoordMat& seeds_coord, std::vector<vec3>& pts) {
  pts.resize(seeds_coord.size());
  for(size_t i = 0; i < pts.size(); ++i) {
    pts[i].x = seeds_coord[i][0];
    pts[i].y = seeds_coord[i][1];
    pts[i].z = seeds_coord[i][2];
  }
}

void get_RDT(RestrictedVoronoiDiagram& RVD,
    dzw::common::VertexCoordMat& rdt_vertex,
    dzw::common::VertexOfFaceMat& rdt_face)
{
  rdt_vertex.resize(RVD.delaunay()->nb_vertices());
  for(unsigned int i = 0; i < RVD.delaunay()->nb_vertices(); ++i) {
    rdt_vertex[i].resize(3, 1);
    for(size_t j = 0; j < rdt_vertex[i].size(1); ++j)
      rdt_vertex[i][j] = RVD.delaunay()->vertex(i)[j];
  }
  RVD.for_each_primal_triangle(get_primal_triangle(rdt_face)) ;
}

void adjust_RVD(const dzw::common::VertexCoordMat& vertices,
    const dzw::common::VertexOfFaceMat& face_vec,
    const std::vector<int>& rvd_face_flag,
    dzw::common::VertexCoordMat& rvd_vertex,
    dzw::common::VertexOfFaceMat& rvd_face,
    std::vector<std::vector<int> >& rvd_classify)
{
  std::vector<int> new_index(vertices.size(), -1);
  std::vector<bool> vertices_flag(vertices_flag.size(), true);
  int now_index = 0;
  new_index[0] = now_index++;
  for(size_t i = 1; i < vertices.size(); ++i) {
    for(size_t j = 0; j < i; ++j) {
      if(vertices_flag[j] && norm(vertices[j] - vertices[i]) < dzw::common::nearly_zero) {
        vertices_flag[i] = false;
        new_index[i] = j;
      }
    }
    if(vertices_flag[i])
      new_index[i] = now_index++;
  }
  rvd_vertex.resize(now_index);
  now_index = 0;
  for(size_t i = 0; i < vertices.size(); ++i) {
    if(vertices_flag[i]) {
      rvd_vertex[now_index++] = vertices[i];
    }
  }
  rvd_face.resize(face_vec.size());
  for(size_t i = 0; i < face_vec.size(); ++i) {
    rvd_face[i].resize(face_vec[i].size());
    for(size_t j = 0; j < face_vec[i].size(); ++j) {
      rvd_face[i][j] = new_index[face_vec[i][j]];
    }
  }
  int max_flag = -1;
  for(size_t i = 0; i < rvd_face_flag.size(); ++i) {
    if(max_flag < rvd_face_flag[i])
      max_flag = rvd_face_flag[i];
  }
  rvd_classify.resize(max_flag - 1);
  for(size_t i = 0; i < rvd_face_flag.size(); ++i) {
    rvd_classify[rvd_face_flag[i] - 1].push_back(i);
  }
}

void get_RVD(RestrictedVoronoiDiagram& RVD,
    dzw::common::VertexCoordMat& rvd_vertex,
    dzw::common::VertexOfFaceMat& rvd_face,
    std::vector<std::vector<int> >& rvd_classify)
{
  std::vector<int> rvd_face_flag;
  dzw::common::VertexCoordMat vertices;
  dzw::common::VertexOfFaceMat face_vec;
  bool sym_backup = RVD.symbolic() ;
  RVD.set_symbolic(true) ;
  RVD.for_each_facet(get_RVD_facets(vertices, face_vec, rvd_face_flag)) ;
  RVD.set_symbolic(sym_backup) ;
}

int get_rvd_mesh(const std::string& file_name,
    const dzw::common::VertexCoordMat& seeds_coord,
    dzw::common::VertexCoordMat& rvd_vertex,
    dzw::common::VertexOfFaceMat& rvd_face,
    std::vector<std::vector<int> >& rvd_classify,
    dzw::common::VertexCoordMat& rdt_vertex,
    dzw::common::VertexOfFaceMat& rdt_face)
{
  Mesh M;
  unsigned int nb_borders = M.load(file_name);
  std::vector<vec3> pts;
  set_pts(seeds_coord, pts);
  Delaunay* delaunay = Delaunay::create("CGAL");
  RestrictedVoronoiDiagram RVD(delaunay, &M);

  delaunay->set_vertices(pts);
  get_RVD(RVD, rvd_vertex, rvd_face, rvd_classify);
  get_RDT(RVD, rdt_vertex, rdt_face);
  delete delaunay;
}

//==================================================================================

/**
 * Used by save_RDT().
 */
class SavePrimalTriangle {
public:
    SavePrimalTriangle(
        std::ofstream& out
    ) : out_(&out) {
    }
    void operator()(unsigned int i, unsigned int j, unsigned int k) const {
        (*out_) << "f " << i+1 << " " << j+1 << " " << k+1 << std::endl ;
    }

private:
    std::ofstream* out_ ;
} ;

/**
 * Given a Restricted Voronoi Diagram, saves the Restricted Delaunay
 * Triangulation to a file in alias|wavefront .obj format.
 */
void save_RDT(
    RestrictedVoronoiDiagram& RVD, const std::string& filename
) {
    std::ofstream out(filename.c_str()) ;
    for(unsigned int i=0; i<RVD.delaunay()->nb_vertices(); i++) {
        out << "v " << RVD.delaunay()->vertex(i) << std::endl ;
    }
    RVD.for_each_primal_triangle(SavePrimalTriangle(out)) ;
    out.close();
}

/**
 * used by save_RVD()
 */
class SaveRVDFacets {
public:
    SaveRVDFacets(
        std::ostream& out
    ) : out_(out), cur_v_(1), cur_f_(1) {
    }
    void operator()(unsigned int iv, Mesh* M) const {
        for(unsigned int f=0; f<M->nb_facets(); f++) {
            for(unsigned int i=M->facet_begin(f); i<M->facet_end(f); i++) {
                const vec3& v = M->vertex(i) ;
                out_ << "v " << v << std::endl ;
            }
            out_ << "f " ;
            for(unsigned int i=M->facet_begin(f); i<M->facet_end(f); i++) {
                out_ << cur_v_ << " ";
                cur_v_++ ;
            }
            out_ << std::endl ;
            cur_f_++ ;
        }
    }
private:
    mutable std::ostream& out_ ;
    mutable unsigned int cur_v_ ;
    mutable unsigned int cur_f_ ;
} ;

/**
 * Saves a Restricted Voronoi Diagram to a file in alias|wavefront .obj format
 * (with Graphite extensions: facet attributes, rename as .eobj to display).
 */
void save_RVD(RestrictedVoronoiDiagram& RVD, const std::string& filename) {
    std::ofstream out(filename.c_str()) ;
    if(!out) {
        std::cerr << "could not open file." << std::endl ;
        return ;
    }
    bool sym_backup = RVD.symbolic() ;
    RVD.set_symbolic(true) ;
    RVD.for_each_facet(SaveRVDFacets(out)) ;
    RVD.set_symbolic(sym_backup) ;
}

int save_rdt_mesh(const std::string& base_file_name,
    const std::string& rvd_file_name,
    const std::string& rdt_file_name,
    const dzw::common::VertexCoordMat& seeds_coord)
{
  Mesh M;
  unsigned int nb_borders = M.load(base_file_name);
  std::vector<vec3> pts;
  set_pts(seeds_coord, pts);
  Delaunay* delaunay = Delaunay::create("CGAL");
  RestrictedVoronoiDiagram RVD(delaunay, &M);

  delaunay->set_vertices(pts);
  save_RVD(RVD, rvd_file_name);
  save_RDT(RVD, rdt_file_name);
  delete delaunay;
}

//================================================================

/**
 * Used by get_combinatorics() in volume mode
 */
class MemorizeIndices {
public:
    MemorizeIndices(
        std::vector<int>& I_in,
        std::vector<vec3>& C_in
    ) : I(I_in), C(C_in) {
        I.resize(0) ;
        C.resize(0) ;
    }

    void operator() (
        unsigned int i,
        int j,
        const VertexEdge& v1,
        const VertexEdge& v2,
        const VertexEdge& v3
    ) const {
        I.push_back(i) ;
        I.push_back(v1.sym[2]) ;
        I.push_back(v1.sym[1]) ;
        I.push_back(v1.sym[0]) ;
        I.push_back(v2.sym[2]) ;
        I.push_back(v2.sym[1]) ;
        I.push_back(v2.sym[0]) ;
        I.push_back(v3.sym[2]) ;
        I.push_back(v3.sym[1]) ;
        I.push_back(v3.sym[0]) ;
        C.push_back(v1) ;
        C.push_back(v2) ;
        C.push_back(v3) ;
    }
private:
    mutable std::vector<int>& I ;
    mutable std::vector<vec3>& C ;
} ;

/**
 * Used by get_combinatorics() in surface mode
 */
class MemorizeIndicesAndFacets{
public:
    MemorizeIndicesAndFacets(
        const RestrictedVoronoiDiagram& RVD_in,
        std::vector<int>& I_in,
        std::vector<vec3>& C_in,
        std::vector<int>& F_in
    ) : RVD(RVD_in), I(I_in), C(C_in), F(F_in) {
        I.resize(0) ;
        C.resize(0) ;
        F.resize(0) ;
    }

    void operator() (
        unsigned int i,
        const VertexEdge& v1,
        const VertexEdge& v2,
        const VertexEdge& v3
    ) const {
        I.push_back(i) ;
        I.push_back(v1.sym[2]) ;
        I.push_back(v1.sym[1]) ;
        I.push_back(v1.sym[0]) ;
        I.push_back(v2.sym[2]) ;
        I.push_back(v2.sym[1]) ;
        I.push_back(v2.sym[0]) ;
        I.push_back(v3.sym[2]) ;
        I.push_back(v3.sym[1]) ;
        I.push_back(v3.sym[0]) ;
        F.push_back(RVD.current_facet()) ;
        C.push_back(v1) ;
        C.push_back(v2) ;
        C.push_back(v3) ;
    }
private:
    const RestrictedVoronoiDiagram& RVD ;
    mutable std::vector<int>& I ;
    mutable std::vector<vec3>& C ;
    mutable std::vector<int>& F ;
} ;


/**
 * Gets the combinatorics of the integration simplices,
 * i.e. 10 integers per integration simplex.
 * (see Section 3.1 in the paper)
 * Returns also the array of C vertices (three per integration simplex).
 * Since they are easy to get during the combinatorial phase, they are
 * computed here and kept for the algebraic phase.
 *
 * In 2D mode (volume = false), returns also the array F.
 *   F[i] indicates the facet that contains the i-th integration simplex.
 *
 */
void get_combinatorics(
    Mesh* M, const std::vector<vec3>& pts,
    std::vector<int>& I, std::vector<vec3>& C, std::vector<int>& F
) {
    Delaunay* delaunay = Delaunay::create("CGAL") ;std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    delaunay->set_vertices(pts) ;std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    RestrictedVoronoiDiagram RVD(delaunay,M) ;std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    RVD.set_symbolic(true) ;std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    RVD.for_each_triangle(MemorizeIndicesAndFacets(RVD,I,C,F)) ;std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    delete delaunay ;
}

/**
 * Computes F_{L_p} and its gradient.
 */
double compute_F_g(Mesh* m, const std::vector<vec3>& pts, std::vector<double>& gradient) {
  std::vector<int> I ;
  std::vector<vec3> C ;
  std::vector<int> F ;
  get_combinatorics(m, pts, I, C, F) ;
  std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  unsigned int nb_integration_simplices = (unsigned int)I.size() / 10 ;
  std::vector<mat3> M(nb_integration_simplices) ;
  for(unsigned int i=0; i<M.size(); i++) {
    M[i].load_identity() ;
    // or replace with anisotropy field
    //   In 2D: use F[i] to retreive the index of the facet that contains
    //      the current integration simplex (and access an array of per-facet anisotropy).
    //   In 3D: use geometric search from the centroid of the current
    //      integration simplex.
  }
  std::vector<plane3> Q(m->nb_facets()) ;
  for(unsigned int i=0; i<m->nb_facets(); i++) {
    Q[i] = m->facet_plane(i) ;
  }
  gradient.resize(pts.size() * 3) ;
  return compute_F_Lp(false, 2, m, I, C, pts, Q, M, gradient) ;
}

void load_pts(const std::string& filename, std::vector<vec3>& pts) {
    pts.clear() ;
    std::ifstream in_stream(filename.c_str()) ;
    if(!in_stream) {
        std::cerr << "Could not open " << filename << std::endl ;
        return ;
    }
    Geex::LineInputStream in(in_stream) ;
    while(!in.eof()) {
        in.get_line() ;
        std::string kw ;
        in >> kw ;
        if(kw == "v") {
            vec3 v ;
            in >> v ;
            pts.push_back(v) ;
        }
    }
}
/*
double get_rvd_f_g(const std::string& file_name,
    const std::string& pts_file)
{
  Mesh M ;
  unsigned int nb_borders = M.load(file_name) ;
  std::vector<vec3> pts ;
  load_pts(pts_file, pts);
  std::vector<double>   gradient;
  return compute_F_g(&M, pts, gradient) ;
}
*/
double get_rvd_f_g(const std::string& file_name,
    const dzw::common::VertexCoordMat& seeds_coord,
    std::vector<double>& gradient)
{
  Mesh M ;
  unsigned int nb_borders = M.load(file_name) ;
  std::vector<vec3> pts ;
  // Convert the format of seed coordinate
  set_pts(seeds_coord, pts);
  return compute_F_g(&M, pts, gradient) ;
}
