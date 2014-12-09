#ifndef HJ_SYMMETRY_HALF_EDGE_H_
#define HJ_SYMMETRY_HALF_EDGE_H_

#include <hjlib/half_edge/half_edge.h>
#include <hjlib/half_edge/container.h>
#include <hjlib/half_edge/operation.h>
#include <vector>

#include "print_half_edge.h"
#include "validate_half_edge.h"

using namespace std;

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

//! @brief return one of the edge from "from" to "to" in "m".
template <typename TMPL, typename CVI>
typename TMPL::edges_t::const_iterator
find_edge(hj::half_edge::mesh_tmpl<TMPL> &m, const CVI &from, const CVI &to)
{
  typedef typename TMPL::edges_t::const_iterator CEI;
  /*
  for(CEI ei = to->edge();!!ei;) {
    assert(ei->vert() == to);//why == to?
    if(ei->oppo()->vert() == from)
      return ei;
    if(ei->next() == CEI()) { // to is an end vertex of an edge
      if(ei->oppo()->vert() == from)
        return ei;
      else
        return CEI();
    }
    ei = ei->next()->oppo();
    if(ei == to->edge())
      break;
      }*/

  //An assertion before find edge
  //first, there is no isolated edge which doesn't belongs to a face
  //second, if an edge belongs to a face, its next and prev are not null
  for (CEI ei = to->edge(); !!ei;)
  {
    //assert(ei->vert());
    if (ei->vert() == from)
    {
      return ei;
    }

    ei = ei->oppo()->next();
    if (ei == to->edge())
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
  for(auto ei = vi->edge();!!ei;) {
    cout << "one-ring edge: " << print_edge(ei, id_key()) << endl;
    if(!ei->face()) { // an in bd
      cout << "add isec: " << endl;
      sec.push_back(ei);
    }
    if(!ei->oppo()->face()) { // an out bd
      cout << "add osec: " << endl;
      sec.push_back(ei->oppo());
    }
    ei = ei->oppo()->next();//MODIFIED prev()-->next()
    if(ei == vi->edge())
      break;
  }
  assert(sec.size()%2 == 0);
  if(!sec.empty() && sec[0]->vert() != vi) {
    std::rotate(sec.begin(), sec.begin()+1, sec.end());
  }
#if 1 // condition for sectors.
  for(size_t i = 0; i < sec.size()/2; ++i) {
    assert(sec[i*2]->vert() == vi);
    assert(sec[i*2+1]->oppo()->vert() == vi);
    //    assert(sec[i*2]->next() == sec[(sec.size()+i*2-1)%sec.size()]);
    //    assert(sec[i*2+1]->prev() == sec[(sec.size()+i*2+2)%sec.size()]);
  }
#endif
}

//! @brief add a new sector from in to out into sec., in and out are
//! outer edge of a sector, but may be not border edges.
template <typename CEI, typename Con>
int merge_sectors(Con &sec, const CEI &in, const CEI &out)
{
  assert(in->vert() == out->oppo()->vert());
  assert(sec.size()%2 == 0);
  size_t i;
  if(!in->face() && !out->face()) { // add a new sector
    sec.push_back(in);
    sec.push_back(out);
  }
  if(!!in->face() && !!out->face()) { // remove a gap
  }
  if(!in->face() && !!out->face()) { // grow a sector
    for(i = 0; i < sec.size()/2; ++i) {
      if(sec[i*2]->oppo() == out) {
        sec[i*2] = in;
      }
    }
    if(i == sec.size()/2)
      return __LINE__;
  }
  if(!!in->face() && !out->face()) { // grow a sector
    for(i = 0; i < sec.size()/2; ++i) {
      if(sec[i*2+1]->oppo() == in) {
        sec[i*2+1] = out;
      }
    }
    if(i == sec.size()/2)
      return __LINE__;
  }
}

template <typename TMPL, typename RND_ITR>
typename TMPL::faces_t::const_iterator
add_face2(hj::half_edge::mesh_tmpl<TMPL> &m, const RND_ITR &vert_loop_beg, std::size_t n)
{
  assert(n > 2); // TODO: may be also good for n=1,2.
  for(size_t i = 0; i < n; ++i) {
    cout << "one-ring-0: " << print_one_ring(vert_loop_beg[i], id_key()) << endl;
  }
  typedef typename TMPL::faces_t::const_iterator CFI;
  CFI fi = m().add(typename TMPL::face_t());

  typedef typename TMPL::edges_t::const_iterator CEI;
  typedef std::vector<CEI> e_con;
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
  // set prev and next
  vector<bool> is_broken(n, false);
  for(std::size_t i = 0; i < n; ++i) {
    if(!edges[i]->oppo()->next())
      is_broken[i] = true;
    if(!edges[i]->oppo()->prev())//??
      is_broken[(i+1)%is_broken.size()] = true;
  }
  for(std::size_t i = 0; i < n; ++i) { // for all node
    if(!is_broken[i]) continue;
    cout << "broken: " << vert_loop_beg[i]->id() << endl;
    assert(edges[i]->oppo()->vert() == vert_loop_beg[i]);
    // get all bd edges at the node vert_loop_beg[i])
    std::vector<CEI> bd; // in, out, in, out, ...
    bd.reserve(6); // rough estimation for efficiency.
    sectors(vert_loop_beg[i], bd);

    // update sectors
    merge_sectors(bd, edges[i]->oppo(), edges[(n+i-1)%n]->oppo());
    // if(!edges[i]->oppo()->face()) { // not cancel
    //   bd.push_back(edges[i]->oppo());
    // }
    // if(!edges[(n+i-1)%n]->oppo()->face())
    //   bd.push_back(edges[(n+i-1)%n]->oppo());

    for(size_t i = 0; i < bd.size()/2; ++i) {
      cout << print_edge(bd[i*2], id_key()) << " " << print_edge(bd[i*2+1], id_key()) << "; ";
    }
    cout << endl;
    // assume ibd[i]->obd[i] is a filled sector.
    // set next and prev for bd edges
    for(size_t j = 0; j < bd.size()/2; ++j) {
      m()[bd[j*2]].next() = bd[(bd.size()+j*2-1)%bd.size()];
      m()[bd[(bd.size()+j*2-1)%bd.size()]].prev() = bd[j*2];
    }
  }
  //  if(g_debug)
  //    exit(0);
  // set edge->face and vert->edge
  for(std::size_t i = 0; i < n; ++i) {
    m()[edges[i]].next() = edges[(i+1)%n];
    m()[edges[i]].prev() = edges[(i+n-1)%n]; // be careful of size_t
    if(!vert_loop_beg[i]->edge()) // first edge for vert vert_loop_beg[i]
      m()[vert_loop_beg[i]].edge() = edges[(n+i-1)%n]; // ei points to i
    m()[edges[i]].face() = fi; // the edge adjacent to this face
    assert(vert_loop_beg[i]->edge()->vert() == vert_loop_beg[i]);
  }
  if(g_debug)
    exit(0);
  for(std::size_t i = 0; i < n; ++i) {
    cout << "one-ring-1: " << print_one_ring(vert_loop_beg[i], id_key()) << endl;
  }
  assert(edges[0]);
  m()[fi].edge() = edges[0];
  return fi;
}

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

#endif
