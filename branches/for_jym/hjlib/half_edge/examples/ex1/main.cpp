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
#include <zjucad/matrix/matrix.h>

using namespace std;

using namespace hj::half_edge;

typedef zjucad::matrix::matrix<double> pt_coord;


struct internal_property
{
  struct vert_t {
    vert_t():id_(v_id++)
    {std::cout << "new vert : " << v_id << std::endl; }
             
    vert_t(size_t id): id_(id){}

    void assign(const vert_t& x)
    {
      id_ = x.id_; coord_ = x.coord_;
    }

    size_t id_;
    pt_coord coord_;
  };
  struct edge_t {
    edge_t():id_(e_id++), split_info({-1,0})
    {std::cout << "new edge : " << e_id << std::endl; }

    edge_t(size_t id): id_(id){}

    void assign(const edge_t& x)
    {
      id_ = x.id_; len_ = x.len_;
    }

    size_t id_;
    double len_;
    struct Split_info
    {
      
      int root;
      size_t level;
    } split_info;
  };
  struct face_t {
    face_t():id_(f_id++)
    {std::cout << "new face : " << f_id << std::endl; }

    face_t(size_t id): id_(id){}

    void assign(const face_t& x)
    {
      id_ = x.id_; normal_ = x.normal_;
    }
    
    size_t id_;
    vector<double> normal_;
  };
  static size_t v_id;
  static size_t e_id;
  static size_t f_id;
};

template <typename TMPL>
void mesh_info(const mesh_tmpl<TMPL> &mt)
{
  typedef TMPL mesh_t;
  const mesh_t &m = mt();
  cout << "============ vert ===========" << endl;
  typedef typename mesh_t::verts_t::const_iterator cvi_t;

  for(cvi_t vi = m.verts().begin(); vi != m.verts().end(); ++vi) {
    //cout << __LINE__ << endl;
    if (!vi->edge())
    {
      break;
    }
    cout << "(" << vi->id_ << ", "
         << vi->edge()->oppo()->vert()->id_  << ")" << endl;
    
    typename mesh_t::edges_t::const_iterator ei = vi->edge();
    for(;!!ei;) {
      assert(ei);
      assert(ei->vert());
      
      cout << "face:";
      if (ei->face()) cout << ei->face()->id_ << ", ";
      else cout << "null, ";

      cout << "edge:";
      cout << ei->oppo()->vert()->id_ << "->" << vi->id_;
      cout << " cur:" << ei->id_
           << " p->" << ei->prev()->id_
           << " n->" << ei->next()->id_ << "), "<< endl;;
      ei = ei->next()->oppo();
      if(ei == vi->edge())
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
size_t set_edge_id(mesh_tmpl<MESH>& mt, size_t &max_id)
{
  typedef MESH mesh_t;
  mesh_t &m = mt();

  typedef typename mesh_t::edges_t::const_iterator edge_itr_t;
  typedef typename mesh_t::verts_t::const_iterator vert_itr_t;
  typedef typename mesh_t::faces_t::const_iterator face_itr_t;

  for (edge_itr_t it = m.edges().begin(); it != m.edges().end(); ++it) {
    if (it->id_ != (size_t)-1)
      max_id = std::max(it->id_, max_id);
  }

  for (edge_itr_t it = m.edges().begin(); it != m.edges().end(); ++it) {
    if (it->id_ == (size_t)-1)
      m[it].id_ = max_id; ++max_id;
  }
  return max_id;
}

template <template <typename> class ICon>
int build_mesh_with_boundary(void)
{
  typedef half_edge_mesh_t<ICon, internal_property> mesh_t;
  typedef typename mesh_t::edges_t::const_iterator edge_itr_t;
  typedef typename mesh_t::verts_t::const_iterator vert_itr_t;
  typedef typename mesh_t::faces_t::const_iterator face_itr_t;

  cout << "\n===== Test mesh with Boundary and arbitrary polygons ====" << endl;
  mesh_t m;

  /*
  std::cout << "-- add vertices" << std::endl;
  typename mesh_t::verts_t::const_iterator v_itrs[4];
  for(size_t i = 0; i < 4; ++i) {
    v_itrs[i] = m.add(typename mesh_t::vert_t(i));
  }

 
  size_t table[] = {
    0, 1, 2,
    2, 1, 3,
    0, 2, 3
  };

  size_t st[4] = {0, 3, 6, 9};

  typename mesh_t::faces_t::const_iterator f_itrs[3];
  
  size_t f = 0;
  size_t e_cnt = 0;
  std::cout << "-- add faces" << std::endl;
  for (size_t f=0; f < 3; ++f) {
    std::cout << "   add face: " << f << std::endl;
    vector<typename mesh_t::verts_t::const_iterator> vert_loop;
    for(size_t e = 0; e < st[f+1]-st[f]; ++e) // the four vertices
    {
      vert_loop.push_back(v_itrs[table[st[f]+e]]);
    }
    
    f_itrs[f] = add_face(m, vert_loop.begin(), st[f+1]-st[f]);

    m[f_itrs[f]].id_ = f;

    for(typename mesh_t::edges_t::const_iterator i =
            f_itrs[f]->edge();;) {
      assert(i);
      m[i].id_ = e_cnt; ++e_cnt;
      i = i->next();
      if(i == f_itrs[f]->edge()) break;
    }
  }
  /*/// add vertices
  std::cout << "-- add vertices" << std::endl;
  typename mesh_t::verts_t::const_iterator v_itrs[15];
  for(size_t i = 0; i < 15; ++i) {
    v_itrs[i] = m.add(typename mesh_t::vert_t());
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
    5, 12, 6,
    4, 5, 6,
    6, 12, 11
  };
  size_t st[13] = {0, 3, 6, 10, 15, 18, 21, 24, 27, 30, 33, 36, 39};

  typename mesh_t::faces_t::const_iterator f_itrs[12];
  typedef typename mesh_t::edges_t::const_iterator CEI;
  typedef typename mesh_t::verts_t::const_iterator CVI;
  typedef typename mesh_t::faces_t::const_iterator CFI;
  
  size_t f = 0;
  size_t e_cnt = 0;
  std::cout << "-- add faces" << std::endl;
  for (size_t f=0; f < 12; ++f) {
    std::cout << "   add face: " << f << std::endl;
    vector<typename mesh_t::verts_t::const_iterator> vert_loop;
    for(size_t e = 0; e < st[f+1]-st[f]; ++e) // the four vertices
    {
      vert_loop.push_back(v_itrs[table[st[f]+e]]);
    }
    f_itrs[f] = add_face2(m, vert_loop.begin(), st[f+1]-st[f]);
  }
  //*/
  mesh_info(m);

  std::cout << "before split edge" << std::endl;
  /*
  CEI ei=m().edges().begin();
  for (; !ei->face(); ++ei);
  std::cout << "edge: " << ei->id_ << std::endl;
  std::cout << "face: " << ei->face()->id_ << std::endl;
  del(m, ei->face());
  //*/
  //*
  CEI ei = m().edges().begin();
  vector<edge_itr_t> edges;
  edges.push_back(ei);
  edges.push_back(ei->next());
  edges.push_back(ei->prev());
  split_edges(m, edges.begin(), 3);
  mesh_info(m);

  //*/
  return 0;
}

size_t internal_property::v_id = 0;
size_t internal_property::e_id = 0;
size_t internal_property::f_id = 0;
int main(int argc, char *argv[])
{
  std::cout << "\n################# TEST std_list  ######################" << std::endl;
  build_mesh_with_boundary<std_list>();

  std::cout << "\n################# TEST std_vector  ####################" << std::endl;
  //build_mesh_with_boundary<std_vector>();
  
  return 0;
}
