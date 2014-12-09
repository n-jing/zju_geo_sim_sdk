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

using namespace std;

using namespace hj::half_edge;

struct internal_property
{
  struct vert_t {
    vert_t():id_(-1){}
    vert_t(size_t id): id_(id){}

    void assign(const vert_t& x)
    {
      id_ = x.id_; pos_ = x.pos_;
    }

    size_t id_;
    vector<double> pos_;
  };
  struct edge_t {
    edge_t():id_(-1){}
    edge_t(size_t id): id_(id){}

    void assign(const edge_t& x)
    {
      id_ = x.id_; len_ = x.len_;
    }

    size_t id_;
    double len_;
  };
  struct face_t {
    face_t():id_(-1){}

    face_t(size_t id): id_(id){}

    void assign(const face_t& x)
    {
      id_ = x.id_; normal_ = x.normal_;
    }
    
    size_t id_;
    vector<double> normal_;
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
    cout << ci->id_ << ": ";
    for(typename mesh_t::edges_t::const_iterator ei = ci->edge();;) {
      assert(ei);
      assert(ei->vert());

      if (ei->face()) cout << '(' << ei->face()->id_ << ", ";
      else cout << "(null, ";
      cout << ei->oppo()->vert()->id_ << "), ";
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
      cout << ei->vert()->id_ << (!ei->oppo()?'t':'f') << " ";
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

template <template <typename> class ICon>
int build_mesh_with_boundary(void)
{
  cout << "\n===== Test mesh with split edge function ====" << endl;
  typedef half_edge_mesh_t<ICon, internal_property> mesh_t;
  mesh_t m;

  // add vertices
  std::cout << "-- add vertices" << std::endl;

  const size_t vs = 15;
  const size_t fs = 12;

  typename mesh_t::verts_t::const_iterator v_itrs[vs];
  for(size_t i = 0; i < vs; ++i) {
    v_itrs[i] = m.add(typename mesh_t::vert_t(i));
  }

  //  add faces: fs faces
  // size_t table[] = {
  //     0, 2, 1,
  //     0, 1, 4,
  //     0, 3, 2,
  //     3, 0, 4,
  //     1, 2, 4,
  //     2, 3, 4,
  // };
  // size_t st[fs+1] = {0, 3, 6, 9, 12, 15, 18};

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

  typename mesh_t::faces_t::const_iterator f_itrs[fs];

  typedef typename mesh_t::edges_t::const_iterator edge_itr_t;
  typedef typename mesh_t::verts_t::const_iterator vert_itr_t;
  typedef typename mesh_t::faces_t::const_iterator face_itr_t;
  
  size_t f = 0;
  size_t e_cnt = 0;
  std::cout << "-- add faces" << std::endl;
  for (size_t f=0; f < fs; ++f) {
    std::cout << "   add face: " << f << std::endl;
    vector<typename mesh_t::verts_t::const_iterator> vert_loop;
    vert_loop.clear();
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
  mesh_info(m);

  // test the split function
  std::cout << "# Split the edge(two vert): ";
  edge_itr_t e = m.edges().begin();

  vert_itr_t new_vert = split_edge(m, e);
  if (new_vert) m()[new_vert].id_ = vs+1;
  ++e;
  vert_itr_t new_vert2 = split_edge(m, e);
  if (new_vert2) m()[new_vert2].id_ = vs+2;
  e = get_edge<mesh_t>(v_itrs[11], v_itrs[14]);
  vert_itr_t new_vert3 = split_edge(m, e);
  if (new_vert3) m()[new_vert3].id_ = vs+3;
  e = get_edge<mesh_t>(v_itrs[7], v_itrs[10]);
  vert_itr_t new_vert4 = split_edge(m, e);
  if (new_vert4) m()[new_vert4].id_ = vs+4;

  std::cout << "# After spliting edge... " << std::endl;
  mesh_info(m);
  
  std::cout << "# After split edge" << std::endl;
  std::cout << "Gain a new vert(ID): " << new_vert->id_ << std::endl;
  mesh_info(m);
  
  return 0;
}

int main(int argc, char *argv[])
{
  std::cout << "\n################# TEST std_list  ######################" << std::endl;
  build_mesh_with_boundary<std_list>();

  std::cout << "\n################# TEST std_vector  ####################" << std::endl;
  build_mesh_with_boundary<std_vector>();
  
  return 0;
}
