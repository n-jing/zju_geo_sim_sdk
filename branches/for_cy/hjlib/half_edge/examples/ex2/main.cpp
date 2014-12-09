#include "../../include/half_edge.h"
#include "../../include/container.h"
#include "../../include/operation.h"

#include <chrono> // c++11

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>

#include <jtflib/mesh/io.h>

using namespace std;
using namespace zjucad::matrix;
using namespace hj::half_edge;

class entity_with_id
{
public:
  virtual ~entity_with_id(){}
  operator size_t() const { return id_; }
  operator size_t&() { return id_; }
  size_t id_;
};

void set_id(entity_with_id &e, size_t id)
{
  e.id_ = id;
}

struct internal_property
{
  struct vert_t : public entity_with_id {
  };
  struct edge_t {
  };
  struct face_t : public entity_with_id {
  };
};

template <typename TMPL, typename ME>
int
half_edge_from_mat(mesh_tmpl<TMPL> &mt, const zjucad::matrix::matrix_expression<ME> &cells, size_t vert_num = -1)
{
  using namespace zjucad::matrix;
  typedef TMPL mesh_t;
  mesh_t &m = mt();
  vert_num = (vert_num == -1)?(zjucad::matrix::max(cells())+1):vert_num;
  matrix<typename mesh_t::verts_t::const_iterator> vis(vert_num);
  for(size_t vi = 0; vi < vert_num; ++vi) {
    vis[vi] = m.add(typename mesh_t::vert_t());
    set_id(m[vis[vi]], vi);
  }
  for(size_t fi = 0; fi < cells().size(2); ++fi) {
    const matrix<typename mesh_t::verts_t::const_iterator> vl
      = vis(cells()(colon(), fi));
    auto f = add_face(m, vl.begin(), vl.size());
    set_id(m[f], fi);
  }
  if (set_opposite_and_boundary_edge(m)) {
    std::cout << "# set opposite and boundary edge: ERROR, can not build mesh." << std::endl;
    return 1;
  }
  return 0;
}

template <typename II, typename IO>
void
copy_i(const II &i_beg, const II &i_end, IO o_beg)
{
  for(II i = i_beg; i != i_end; ++i, ++o_beg) *o_beg = i;
}

class timer
{
public:
  timer(const char *msg):beg_(hrc_.now()), msg_(msg){}
  ~timer() {
    cout << msg_ << ": " << chrono::duration_cast<chrono::milliseconds>(hrc_.now()-beg_).count() << " ms" << endl;
  }
private:
  chrono::high_resolution_clock hrc_;
  chrono::high_resolution_clock::time_point beg_;
  const string msg_;
};

template <template <typename> class ICon>
int performance(const matrix<double> &nods, const matrix<size_t> &tris)
{
  typedef half_edge_mesh_t<ICon, internal_property> mesh_t;
  mesh_t m;
  {
    timer t("half_edge_from_mat");
    if (half_edge_from_mat(m, tris)) return 1;
  }
  {
    cout << "=== vertices in the first face ===" << endl;
    auto f0 = m.faces().begin();
    cout << "face id: " << *f0 << endl;
    for(auto ei = f0->edge();;) {
      auto vi = ei->vert();
      cout << "vert id: " << *vi << " with pos " << trans(nods(colon(), *vi));
      ei = ei->next();
      if(ei == f0->edge())
        break;
    }
  }
  vector<typename mesh_t::faces_t::const_iterator> fis(m.faces().size());
  {
    timer t("gather iterator for random access");
    copy_i(m.faces().begin(), m.faces().end(), fis.begin());
  }
  matrix<size_t> face_to_del = rand<double>(tris.size(2)*0.2, 1)*(tris.size(2)-1);
  {
    timer t("is_valid");
    for(size_t i = 0; i < face_to_del.size()/10; ++i)
      is_valid(m, fis[face_to_del[i]]);
  }
  sort(face_to_del.begin(), face_to_del.end());
  const size_t num = unique(face_to_del.begin(), face_to_del.end()) - face_to_del.begin();
  face_to_del = temp(face_to_del(colon(0, num-1)));
  random_shuffle(face_to_del.begin(), face_to_del.end());
  {
    timer t("random delete face");
    for(size_t i = 0; i < face_to_del.size(); ++i)
      del(m, fis[face_to_del[i]]);
  }
  return 0;
}

int main(int argc, char *argv[])
{
  if(argc < 2) {
    cerr << "Usage: ex2 obj" << endl;
    return __LINE__;
  }

  matrix<double> nods;
  matrix<size_t> tris;

  if(jtf::mesh::load_obj(argv[1], tris, nods)) {
    cerr << "load " << argv[1] << " fail." << endl;
    return __LINE__;
  }

  std::cout << "\n################# TEST std_list  ######################" << std::endl;
  performance<std_list>(nods, tris);

  std::cout << "\n################# TEST std_vector  ####################" << std::endl;  
  performance<std_vector>(nods, tris);
  
  return 0;
}
