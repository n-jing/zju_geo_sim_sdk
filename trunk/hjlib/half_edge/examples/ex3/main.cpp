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
      cout << ei->vert()->id_ << "), ";
      ei = ei->prev()->oppo();
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

  #define vs 5
  #define fs 6

  typename mesh_t::verts_t::const_iterator v_itrs[vs];
  for(size_t i = 0; i < vs; ++i) {
    v_itrs[i] = m.add(typename mesh_t::vert_t(i));

  }

  //  add faces: fs faces
  size_t table[] = {
      0, 2, 1,
      0, 1, 4,
      0, 3, 2,
      3, 0, 4,
      1, 2, 4,
      2, 3, 4,
  };
  size_t st[fs+1] = {0, 3, 6, 9, 12, 15, 18};

  typename mesh_t::faces_t::const_iterator f_itrs[fs];


  // #define vs 7
  // #define fs 2
  // typename mesh_t::verts_t::const_iterator v_itrs[vs];
  // for(size_t i = 0; i < vs; ++i) {
  //   v_itrs[i] = m.add(typename mesh_t::vert_t(i));

  // }

  // size_t table[] = {
  //     0, 1, 2, 3, 4,
  //     1, 0, 5, 6,
  // };
  // size_t st[fs+1] = {0, 5, 9};

  // typename mesh_t::faces_t::const_iterator f_itrs[fs];


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
  for (edge_itr_t it = m.edges().begin(); it != m.edges().end(); ++ it) {
      if (it->oppo()->face() != face_itr_t()){
	  e = it;
	  break;
      }
  }

  std::cout << e->vert()->id_ << "---" << e->oppo()->vert()->id_ << std::endl;
  vert_itr_t new_vert = split_edge(m, e);
  if (new_vert == vert_itr_t()) {
      std::cout << "ERROR IN SPLIT_EDGE FUNCTION" << std::endl;
      return 1;
  }

  std::size_t vid = m.verts().size() - 1;
  for (vert_itr_t it = m.verts().begin(); it != m.verts().end(); ++ it) {
      if (!(it->id_ < vid))
	  m()[it].id_ = vid ++;
  }
  
  std::size_t fid = m.faces().size() - 2;
  for (face_itr_t it = m.faces().begin(); it != m.faces().end(); ++ it) {
      if (!(it->id_ < fid))
	  m()[it].id_ = fid ++;
  }

  std::cout << "# After splie edge" << std::endl;
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

