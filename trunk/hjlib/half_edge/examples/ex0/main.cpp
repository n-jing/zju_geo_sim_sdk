#include "../../include/half_edge.h"
#include "../../include/container.h"
#include "../../include/operation.h"

#include <cassert>
#include <cstddef>
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <string>
#include <boost/any.hpp>

using namespace std;

using namespace hj::half_edge;

struct internal_property
{
  //! @note if the entity has copy outside of the mesh, the reflection
  //! will not work
  template <typename T>
  class reflection
  {
  public:
    reflection(const size_t id):id_(id) {
      cout << "ref add: " << id_ << endl;
    }
    ~reflection() {
      cout << "ref del: " << id_ << endl;
    }
    const size_t id_;
  };
  struct vert_t {

    vert_t(){}
    vert_t(size_t id): ref_(new reflection<vert_t>(id)){}
    void assign(const vert_t& x) {
      ref_.reset(new reflection<vert_t>(x.id()));
    }
    
    size_t id() const { return ref_->id_; }
    virtual ~vert_t(){}
  private:
    std::shared_ptr<reflection<vert_t> > ref_;
  };
  struct edge_t {
    edge_t():id_(-1){}
    edge_t(size_t id): id_(id){}
    virtual ~edge_t(){}

    void assign(const edge_t& x) {
      id_ = x.id_;
    }
    size_t id_;
  };
  struct face_t {
    face_t():id_(-1){}
    face_t(size_t id): id_(id){}
    virtual ~face_t(){}

    void assign(const face_t& x) {
      id_ = x.id_;
    }
    size_t id_;
    template <typename T>
    void set(const string &key, const T &v) {
      dyn_prop_[key] = v;
    }
    template <typename T>
    const T &get(const string &key) const {
      map<string, boost::any>::const_iterator i = dyn_prop_.find(key);
      if(i != dyn_prop_.end())
        return boost::any_cast<const T &>(i->second);
      throw;
    }
  private:
    map<string, boost::any> dyn_prop_;
  };
};

template <typename TMPL>
void mesh_info(const mesh_tmpl<TMPL> &mt)
{
  typedef TMPL mesh_t;
  const mesh_t &m = mt();
  cout << "============ vert ===========" << endl;
  typedef typename mesh_t::verts_t::const_iterator cvi_t;
  for(cvi_t ci = m.verts().begin(); ci != m.verts().end(); ++ci) {
    //cout << __LINE__ << endl;
    cout << ci->id() << ": ";
    for(typename mesh_t::edges_t::const_iterator ei = ci->edge();;) {
      assert(ei);
      assert(ei->vert());
      if (ei->face()) cout << '(' << ei->face()->id_ << ", ";
      else cout << "(null, ";
      cout << ei->oppo()->vert()->id() << "), ";
      ei = ei->next()->oppo();
      if(ei == ci->edge())
        break;
    }
    cout << endl;
  }
  cout << "============ face ===========" << endl;
  typedef typename mesh_t::faces_t::const_iterator cfi_t;
  for(cfi_t ci = m.faces().begin(); ci != m.faces().end(); ++ci) {
    cout << ci->id_ << ": ";
    for(typename mesh_t::edges_t::const_iterator ei = ci->edge();;) {
      cout << ei->vert()->id() << (!ei->oppo()?'t':'f') << " ";
      ei = ei->next();
      if(ei == ci->edge())
        break;
    }
    cout << endl;
  }
}

template<typename MESH>
void set_edge_id(mesh_tmpl<MESH>& mt)
{
  typedef MESH mesh_t;
  mesh_t &m = mt();

  typedef typename mesh_t::edges_t::const_iterator edge_itr_t;
  typedef typename mesh_t::verts_t::const_iterator vert_itr_t;
  typedef typename mesh_t::faces_t::const_iterator face_itr_t;

  size_t max_id = 0;
  for (edge_itr_t it = m.edges().begin(); it != m.edges().end(); ++it) {
    if (it->id_ != (size_t)-1)
      max_id = std::max(it->id_, max_id);
  }

  for (edge_itr_t it = m.edges().begin(); it != m.edges().end(); ++it) {
    if (it->id_ == (size_t)-1)
      m[it].id_ = max_id; ++max_id;
  }
}

template<typename MESH>
void test_adjacent_info_requirement(const mesh_tmpl<MESH>& mt)
{
  typedef MESH mesh_t;
  const mesh_t &m = mt();

  typedef typename mesh_t::edges_t::const_iterator edge_itr_t;
  typedef typename mesh_t::verts_t::const_iterator vert_itr_t;
  typedef typename mesh_t::faces_t::const_iterator face_itr_t;
  
  
  std::cout << "# Test requiring adjacent info:" << std::endl;
  
  vector<vert_itr_t> fv;
  face_adj_verts<mesh_t>(m.faces().begin(), fv);
  std::cout << "-- face " << m.faces().begin()->id_ << " adj verts: ";
  for (size_t i=0; i<fv.size(); ++i) {
    std::cout << " " << fv[i]->id();
  }
  std::cout << std::endl;
  
  vector<face_itr_t> ff;
  face_adj_faces<mesh_t>(m.faces().begin(), ff);
  std::cout << "-- face " << m.faces().begin()->id_ << " adj faces: ";
  for (size_t i=0; i<ff.size(); ++i) {
    std::cout << " " << ff[i]->id_;
  }
  std::cout << std::endl;

  vector<face_itr_t> vf;
  vert_adj_faces<mesh_t>(m.verts().begin(), vf);
  std::cout << "-- vert " << m.verts().begin()->id() << " adj faces: ";
  for (size_t i=0; i<vf.size(); ++i) {
    std::cout << " " << vf[i]->id_;
  }
  std::cout << std::endl;

  vector<vert_itr_t> vv;
  vert_adj_verts<mesh_t>(m.verts().begin(), vv);
  std::cout << "-- vert " << m.verts().begin()->id() << " adj verts: ";
  for (size_t i=0; i<vv.size(); ++i) {
    std::cout << " " << vv[i]->id();
  }
  std::cout << std::endl;

  vector<edge_itr_t> ve;
  vert_adj_out_edges<mesh_t>(m.verts().begin(), ve);
  std::cout << "-- vert " << m.verts().begin()->id() << " adj out edges: ";
  for (size_t i=0; i<ve.size(); ++i) {
    std::cout << " " << ve[i]->id_;
  }
  std::cout << std::endl;
}

template <template <typename> class ICon>
int build_mesh(void)
{
  cout << "\n===== Test mesh without Boundary ====" << endl;  
  typedef half_edge_mesh_t<ICon, internal_property> mesh_t;
  mesh_t m;

  // add vertices
  typename mesh_t::verts_t::const_iterator v_itrs[8];
  for(size_t i = 0; i < 8; ++i) {
    v_itrs[i] = m.add(typename mesh_t::vert_t(i));
  }

  // add faces
  size_t table[6][4] = {
    0, 1, 2, 3,
    7, 6, 5, 4,
    1, 0, 4, 5,
    2, 1, 5, 6,
    3, 2, 6, 7,
    0, 3, 7, 4
  };

  typename mesh_t::faces_t::const_iterator f_itrs[6];
  for(size_t f = 0; f < 6; ++f) {
    vector<typename mesh_t::verts_t::const_iterator> vert_loop(4);
    for(size_t e = 0; e < 4; ++e) // the four vertices
      vert_loop[e] = v_itrs[table[f][e]];
    f_itrs[f] = add_face(m, vert_loop.begin(), 4);
    m[f_itrs[f]].id_ = f;
    m[f_itrs[f]].set("first node", table[f][0]);
    cout <<  "first node: " << m[f_itrs[f]].get<size_t>("first node") << endl;
    size_t j = 0;
    for(typename mesh_t::edges_t::const_iterator i = f_itrs[f]->edge();; ++j) {
      assert(i);
      m[i].id_ = f*4+j;
      i = i->next();
      if(i == f_itrs[f]->edge())
        break;
    }
  }

  // glue faces
  // vector<typename mesh_t::edges_t::const_iterator> border_edges;
  // find_border(m.edges(), border_edges);
  // glue_2_side_holes(m, border_edges, border_edges.size());
  
  std::cout << "# build halfedge." << std::endl;
  if (set_opposite_and_boundary_edge(m)) {
    cout << "# ERROR: build halfedge." << endl;
    return 1;
  }

  mesh_info(m);

  std::cout << "# delete some faces" << std::endl;
  del(m, f_itrs[0]);
  mesh_info(m);

  std::cout << "\n# Test Copy Mesh" << std::endl;
  mesh_t mb;
  copy(m, mb);
  assert(is_valid(mb));
  mesh_info(mb);

  return 0;
}

template <template <typename> class ICon>
int build_mesh_with_boundary(void)
{
  cout << "\n===== Test mesh with Boundary and arbitrary polygons ====" << endl;
  typedef half_edge_mesh_t<ICon, internal_property> mesh_t;
  mesh_t m;

  // add vertices
  std::cout << "-- add vertices" << std::endl;
  typename mesh_t::verts_t::const_iterator v_itrs[15];
  for(size_t i = 0; i < 15; ++i) {
    v_itrs[i] = m.add(typename mesh_t::vert_t(i));
  }

  // add faces: 12 faces
  size_t table[] = {
    0, 1, 2,
    2, 4, 3,
    3, 4, 6, 7,
    3, 7, 10, 9, 8,
    7, 6, 11,
    7, 11, 10,
    10, 11, 14,
    11, 13, 14,
    11, 12, 13,
    6, 12, 11,
    5, 12, 6,
    4, 5, 6
  };
  size_t st[13] = {0, 3, 6, 10, 15, 18, 21, 24, 27, 30, 33, 36, 39};

  typename mesh_t::faces_t::const_iterator f_itrs[12];
  typedef typename mesh_t::edges_t::const_iterator edge_itr_t;
  typedef typename mesh_t::verts_t::const_iterator vert_itr_t;
  typedef typename mesh_t::faces_t::const_iterator face_itr_t;
  
  size_t f = 0;
  size_t e_cnt = 0;
  std::cout << "-- add faces" << std::endl;
  for (size_t f=0; f < 12; ++f) {
    std::cout << "   add face: " << f << std::endl;
    vector<typename mesh_t::verts_t::const_iterator> vert_loop;
    for(size_t e = 0; e < st[f+1]-st[f]; ++e) // the four vertices
      vert_loop.push_back(v_itrs[table[st[f]+e]]);
    f_itrs[f] = add_face(m, vert_loop.begin(), st[f+1]-st[f]);
    m[f_itrs[f]].id_ = f;
    for(typename mesh_t::edges_t::const_iterator i = f_itrs[f]->edge();;) {
      assert(i);
      m[i].id_ = e_cnt; ++e_cnt;
      i = i->next();
      if(i == f_itrs[f]->edge()) break;
    }
  }

  std::cout << "# build halfedge." << std::endl;
  if (set_opposite_and_boundary_edge(m)) {
    cout << "# ERROR: build halfedge." << endl;
    return 1;
  }

  set_edge_id(m);
  
  cout << "# before delete operation: " << std::endl;
  mesh_info(m);
  test_adjacent_info_requirement(m);
  std::cout << "# delete some faces: 0, 3, 5" << std::endl;
  del(m, f_itrs[0]);
  del(m, f_itrs[3]);
  del(m, f_itrs[5]);
  mesh_info(m);
  test_adjacent_info_requirement(m);
  std::cout << "# delete some verts: the first vert" << std::endl;
  del(m, m.verts().begin());
  del(m, m.verts().begin());
  mesh_info(m);
  test_adjacent_info_requirement(m);
  std::cout << "# delete some edges: the first edge" << std::endl;
  edge_itr_t e = m.edges().begin();
  std::cout << "-- edge vert: " << e->vert()->id() << std::endl;
  std::cout << "-- edge from: " << e->oppo()->vert()->id() << std::endl;
  edge_itr_t v2e = get_edge<mesh_t>(e->vert(), e->oppo()->vert());
  assert((!v2e)==false);
  std::cout << "-- edge of two vert (" << e->vert()->id() << ","
            << e->oppo()->vert()->id() << ") : "
            << v2e->id_ << std::endl;
  del(m, e);
  std::cout << "-- end of del." << std::endl;
  mesh_info(m);
  cout << (is_valid(m)?"valid":"invalid") << endl;

  std::cout << "# Test Valence:" << std::endl;
  std::cout << "-- valence of vert " << m.verts().begin()->id()
            << ": " << valence<mesh_t>(m.verts().begin()) << std::endl;
  std::cout << "-- valence of face " << m.faces().begin()->id_
            << ": " << valence<mesh_t>(m.faces().begin()) << std::endl;

  test_adjacent_info_requirement(m);
  
  std::cout << "\n# Test Copy Mesh" << std::endl;
  mesh_t mb;
  copy(m, mb);
  assert(is_valid(mb));
  mesh_info(mb);
  
  return 0;
}

int main(int argc, char *argv[])
{
  std::cout << "\n################# TEST std_list  ######################" << std::endl;
  build_mesh<std_list>();
  build_mesh_with_boundary<std_list>();

  std::cout << "\n################# TEST std_vector  ####################" << std::endl;  
  build_mesh<std_vector>();
  build_mesh_with_boundary<std_vector>();
  
  return 0;
}
