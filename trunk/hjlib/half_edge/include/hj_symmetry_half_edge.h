#ifndef HJ_SYMMETRY_HALF_EDGE_H_
#define HJ_SYMMETRY_HALF_EDGE_H_

#include <hjlib/half_edge/half_edge.h>
#include <hjlib/half_edge/container.h>
#include <hjlib/half_edge/operation.h>
#include <vector>

//#include "print_half_edge.h"
//#include "validate_half_edge.h"

using namespace std;
/*
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

bool g_debug = false;
//*/
template <typename TMPL, typename CVI>
void
adjust_vert_edge(hj::half_edge::mesh_tmpl<TMPL> &m, const CVI & vert)
{
  //return for condition already holds
  if (!vert->edge()->face())
  {
    cout << "face is null" << endl;
    return;
  }
  
  typedef typename TMPL::edges_t::const_iterator CEI;
  CEI ei = vert->edge();

  while (!!ei)
  {
    if (!(ei->face()) )
    {
      m()[vert].edge() = ei;
    }
    ei = ei->oppo()->next();
    
    if (ei == vert->edge())
    {
      break;
    }
  }
}

//! @brief return one of the edge from "from" to "to" in "m".
template <typename TMPL, typename CVI>
typename TMPL::edges_t::const_iterator
find_edge(hj::half_edge::mesh_tmpl<TMPL> &m, const CVI &from, const CVI &to)
{
  typedef typename TMPL::edges_t::const_iterator CEI;

  for (CEI ei = from->edge(); !!ei;)
  {
    if (ei->vert() == to)
    {
      return ei;
    }

    ei = ei->oppo()->next();
    if (ei == from->edge())
    {
      break;
    }
  }
  return CEI();
}

class id_key
{
 public:
  template <typename CI>
  size_t operator()(const CI &ci) const { return ci->id(); }
};

// @brief all edges in sec are border edge, and the region between
// sec[i*2] and sec[i*2+1] are connected by faces.
template <typename CVI, typename Con>
void sectors(const CVI &vi, Con &sec)
{
  for(auto ei=vi->edge(); ;)
  {
    //    cout << "one-ring edge: " << print_edge(ei, id_key()) << endl;

    if(!ei->face()) {
      cout << "add isec: " << endl;
      //add in and out edges of a sector sequencely
      sec.push_back(ei->prev());
      sec.push_back(ei);
    }

    ei = ei->oppo()->next();
    if(ei == vi->edge())
      break;
  }
}


//! @brief add a new sector from in to out into sec., in and out are
//! outer edge of a sector, but may be not border edges.
template <typename TMPL, typename CEI, typename Con>
int add_face_into_sectors(hj::half_edge::mesh_tmpl<TMPL> &m, Con &sec, const CEI &in, const CEI &out)
{
  assert(in->vert() == out->oppo()->vert());
  assert(sec.size());
  assert(sec.size()%2 == 0);

  bool is_added = false;
  for (size_t i=0; i<sec.size(); ++i)
  {
    //find existing edges
    if (in == sec[i])//only happened in even numbers
    {
      if (out != sec[i+1])
      {
        m()[out->oppo()].next() = sec[i+1];
        m()[sec[i+1]].prev() = out->oppo();
      }
      is_added = true;
    }

    if (out == sec[i])//only happened in odds numbers
    {
      if (in != sec[i-1])
      {
        m()[in->oppo()].prev() = sec[i-1];
        m()[sec[i-1]].next() = in->oppo();
      }
      is_added = true;
    }
  }

  //add the new face to the first sector
  //IMPROVE: this should provides an interface for geometry information
  if (!is_added)
  {
    m()[sec[0]].next() = in->oppo();
    m()[in->oppo()].prev() = sec[0];
    m()[sec[1]].prev() = out->oppo();
    m()[out->oppo()].next() = sec[1];
  }

  adjust_vert_edge(m, in->vert());
}

template <typename TMPL, typename RND_ITR>
typename TMPL::faces_t::const_iterator
add_face2(hj::half_edge::mesh_tmpl<TMPL> &m, const RND_ITR &vert_loop_beg, std::size_t n)
{
  //add vertices
  vector<bool> is_v_exist(n, true);//is the vertex already exists
  for (size_t i=0; i<n; ++i)
  {
    is_v_exist[i] =!!(vert_loop_beg[i]->edge());
  }
  
  assert(n > 2); // TODO: may be also good for n=1,2.

  typedef typename TMPL::faces_t::const_iterator CFI;
  CFI fi = m().add(typename TMPL::face_t());
  typedef typename TMPL::edges_t::const_iterator CEI;
  typedef vector<CEI> e_con;
  e_con edges(n); // edges[i]: v_i -> v_{i+1}

  // add edges, set vert and oppo
  for(size_t i = 0; i < n; ++i) {
    // find edge v_i to v_{i+1}
    CEI ei = find_edge(m(), vert_loop_beg[i], vert_loop_beg[(i+1)%n]);
    if(!ei) { // add the pair of half_edge
      ei = m().add(typename TMPL::edge_t());
      m()[ei].vert() = vert_loop_beg[(i+1)%n];
      CEI eio = m().add(typename TMPL::edge_t());
      m()[eio].vert() = vert_loop_beg[i];
      m()[ei].oppo() = eio;
      m()[eio].oppo() = ei;
    }
    edges[i] = ei;
    assert(!ei->face()); // assert for non-manifold edge
  }

  //modify prev and next for each vertex of the new face
  for (size_t i=0; i<n; ++i)
  {
    //vertex already exist means that their are other edges need to be
    //consider around this vertex
    if (is_v_exist[i])
    {
      vector<CEI> bd;
      bd.reserve(6);
      //note: even after the new edges are added, sector function still
      //can not find these edges because their prev and next are abscent
      sectors(vert_loop_beg[i], bd);
      add_face_into_sectors(m, bd, edges[(i+n-1)%n], edges[i]);
    }
    else
    {
      //new vertex with two new edges, simple topology relation
      m()[edges[(i+n-1)%n]->oppo()].prev() = edges[i]->oppo();
      m()[edges[i]->oppo()].next() = edges[(i+n-1)%n]->oppo();
    }
  }

  // set edge->face and vert->edge
  for(size_t i = 0; i < n; ++i) {
    m()[edges[i]].next() = edges[(i+1)%n];
    m()[edges[i]].prev() = edges[(i+n-1)%n]; // be careful of size_t
    if(!vert_loop_beg[i]->edge()) // first edge for vert vert_loop_beg[i]
      m()[vert_loop_beg[i]].edge() = edges[i]; // ei points to i
    m()[edges[i]].face() = fi; // the edge adjacent to this face
    assert(vert_loop_beg[i]->edge()->oppo()->vert() == vert_loop_beg[i]);
  }

  //adjust edge of every vertex to guarantee boundary vertex have bd edge
  for (size_t j=0; j<n; ++j)
  {
    adjust_vert_edge(m, edges[j]->vert());
  }
  
  assert(edges[0]);
  m()[fi].edge() = edges[0];
  return fi;
}

/*
template <typename TMPL, typename ME>
int
half_edge_from_mat(hj::half_edge::mesh_tmpl<TMPL> &mt, const zjucad::matrix::matrix_expression<ME> &cells, size_t vert_num = -1)
{
  using namespace zjucad::matrix;
  typedef TMPL mesh_t;
  typedef typename mesh_t::verts_t::const_iterator CVI;
  typedef typename mesh_t::edges_t::const_iterator CEI;
  typedef typename mesh_t::faces_t::const_iterator CFI;
  mesh_t &m = mt();
  vert_num = (vert_num == -1)?(zjucad::matrix::max(cells())+1):vert_num;
  matrix<CVI> vis(vert_num);
  for(size_t vi = 0; vi < vert_num; ++vi) {
    vis[vi] = m.add(typename mesh_t::vert_t());
    m[vis[vi]].id() =vi;
  }
  for(size_t fi = 0; fi < cells().size(2); ++fi) {
    const matrix<CVI> vl
        = vis(cells()(colon(), fi));
    std::cout << "add face: " << fi << " " << trans(cells()(colon(), fi)) << std::endl;
    //    g_debug = (fi == 3);
    //    if(fi == 3)
    //      exit(0);
    auto f = add_face2(m, vl.begin(), vl.size());
    m[f].id() = fi;

    //    std::cout << m.verts().size() << " " << m.edges().size() << " " << m.faces().size() << std::endl;
    {
      const CEI ei = find_wrong_edge(m.edges());
      if(!!ei) {
        cout << "error at edge: " << print_edge(ei, id_key()) << endl;
        cout << print_one_ring(ei->vert(), id_key()) << endl
             << print_one_ring(ei->oppo()->vert(), id_key()) << endl;
      }
      pair<CVI, int> err = find_wrong_one_ring(m);
      if(!!err.first)
        cout << "error at valence: " << err.first->id() << " " << err.second << endl;
      const CFI fi = find_wrong_face(m.faces());
      if(!!fi) {
        cout << "error at face: " << fi->id() << endl;
      }
    }
  }
  // if(set_opposite_and_boundary_edge(m)) {
  //   std::cout << "# set opposite and boundary edge: ERROR, can not build mesh." << std::endl;
  //   return 1;
  // }
  return 0;
}
//*/
#endif
