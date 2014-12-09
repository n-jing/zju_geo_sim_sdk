//! @file operation.hh
//! @brief implementation of the basic operations on half_edge,
//! included into operation.h
//! @author Jin Huang, Xianzhong Fang
//! @date 201406-

template <typename TMPL>
bool
is_valid(const cont_tmpl<TMPL> &c, const typename TMPL::const_iterator &i)
{
  typedef TMPL cont_t;
  const cont_t &cont = c();
  for(typename cont_t::const_iterator j = cont.begin(); j != cont.end(); ++j)
    if(i == j)
      return true;
  return false;
}

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m, const typename TMPL::verts_t::const_iterator &vi)
{
  return is_valid(m().verts(), vi) &&
    is_valid(m().edges(), vi->edge());
}

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &ei)
{
  return is_valid(m().edges(), ei) &&
    is_valid(m().edges(), ei->prev()) &&
    is_valid(m().edges(), ei->next()) &&
    is_valid(m().edges(), ei->oppo()) && // the convension of null face for boundary
    is_valid(m().verts(), ei->vert()) &&
    (ei->face()?is_valid(m().faces(), ei->face()):true);
}

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m, const typename TMPL::faces_t::const_iterator &fi)
{
  return is_valid(m().faces(), fi) &&
    is_valid(m().edges(), fi->edge());
}

//! @note may be insert iterator will be better.
template <typename TMPL, typename VC, typename EC, typename FC>
void
find_invalid(const mesh_tmpl<TMPL> &m, VC &vc, EC &ec, FC &fc)
{
  for(auto i = m().verts().begin(); i != m().verts().end(); ++i)
    if(!is_valid(m, i))
      vc.push_back(i);
  for(auto i = m().edges().begin(); i != m().edges().end(); ++i)
    if(!is_valid(m, i))
      ec.push_back(i);
  for(auto i = m().faces().begin(); i != m().faces().end(); ++i)
    if(!is_valid(m, i))
      fc.push_back(i);
}

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m)
{
  std::vector<typename TMPL::verts_t::const_iterator> vc;
  std::vector<typename TMPL::edges_t::const_iterator> ec;
  std::vector<typename TMPL::faces_t::const_iterator> fc;
  find_invalid(m, vc, ec, fc);
  return (vc.size() + ec.size() + fc.size() == 0);
}

//! @brief add a polygon face to mesh without setting opposite edge
template <typename TMPL, typename ITR>
typename TMPL::faces_t::const_iterator
add_face(mesh_tmpl<TMPL> &m, const ITR &vert_loop_beg, std::size_t n)
{
  typename TMPL::faces_t::const_iterator fi = m().add(typename TMPL::face_t());
  typedef std::vector<typename TMPL::edges_t::const_iterator> e_con;
  e_con edges;
  edges.reserve(n);
  ITR vi = vert_loop_beg;
  for(; edges.size() < n; ++vi) {
    typename TMPL::edges_t::const_iterator ei = m().add(typename TMPL::edge_t());
    assert(ei);
    edges.push_back(ei);
    m()[ei].vert() = *vi; // the edge points to this vertex
    m()[ei].face() = fi;  // the edge adjacent to this face
    m()[*vi].edge() = ei; // ei points to vi
  }
  vi = vert_loop_beg;
  for(std::size_t i = 0; i < edges.size(); ++i, ++vi) {
    m()[edges[i]].next() = edges[(i+1)%edges.size()];
    m()[edges[i]].prev() = edges[(i+edges.size()-1)%edges.size()]; // be careful of size_t
  }
  assert(edges[0]);
  m()[fi].edge() = edges[0];
  return fi;
}

//! @brief NOTICE append to border_edges.
template <typename TMPL, typename C>
void
find_border(const cont_tmpl<TMPL> &edges, C &border_edges)
{
  typename C::value_type ci;
  for(ci = edges().begin(); ci != edges().end(); ++ci) {
    if(!(ci->face()))
      border_edges.push_back(ci);
  }
}

//! @brief a boundary vert is a vert with at least one adjacent boundary edge
template <typename TMPL>
bool is_boundary(mesh_tmpl<TMPL>& m, const typename TMPL::verts_t::const_iterator& v)
{
  return !(v->edge()->face());
}
//! @brief a boundary edge is an edge with null face
template <typename TMPL>
bool is_boundary(mesh_tmpl<TMPL>& m, const typename TMPL::edges_t::const_iterator& e)
{
  return !(e->face());
}
//! @brief a boundary face is a face with at least one edge, which has a opposite boundary edge
template <typename TMPL>
bool is_boundary(mesh_tmpl<TMPL>& m, const typename TMPL::faces_t::const_iterator& f)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = f->edge();
  do {
    if (!e->oppo()->face()) return true;
    e = e->next();
  } while (e != f->edge());

  return true;
}

//! @brief an isolated vert is a vert that it doesn't have adjacent edge.
template <typename TMPL>
    bool is_isolated(const typename TMPL::verts_t::const_iterator& v)
{
  return !v->edge();
}
//! @brief an isolated edge is a edge that it and its oppo edge doesn't have adjacent face.
template <typename TMPL>
    bool is_isolated(const typename TMPL::edges_t::const_iterator& e)
{
  return (!e->face() && !e->oppo()->face());
}
//! @brief an isolated face is a face that it doesn't have adjacent face.
template <typename TMPL>
    bool is_isolated(const typename TMPL::faces_t::const_iterator& f)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = f->edge();
  do {
    if (e->oppo()->face()) return false;
    e = e->next();
  } while (e != f->edge());
  return true;
}

//! @brief delete the vertex and the related edges and faces
template <typename TMPL>
int del(mesh_tmpl<TMPL> &m, const typename TMPL::verts_t::const_iterator& v)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;

  vert_itr_t vv = v; // copy because v is a reference
  
  std::vector<face_itr_t> vf;
  vert_adj_faces<TMPL>(vv, vf);
  for (auto it = vf.begin(); it != vf.end(); ++it) {
    if (del(m, *it)) return 1;
  }

  std::vector<edge_itr_t> ve;
  vert_adj_out_edges<TMPL>(vv, ve);
  for (auto it = ve.begin(); it != ve.end(); ++it) {
    if (del(m, *it)) return 1;
  }

  m().del(vv);

  return 0;
}

//! @brief delete the edge through the related faces
template <typename TMPL>
int del(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator& e)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;

  edge_itr_t ee = e; // copy, e is a reference, the next operation may modify it
  auto ee_op = ee->oppo();
  if (ee->face()) {
    if (del(m, ee->face())) return 1;
  }
  if (ee_op->face()) { // can't use e->oppo()->face() directly, because the e may be deleted
    if (del(m, ee_op->face())) return 1;
  }

  // adjust the edge that vert point to
  vert_itr_t v = ee->vert();
  if (v->edge() == ee) {
    m()[v].edge() = ee->oppo()->prev();
    if (v->edge() == ee)
    {
      m()[v].edge() = edge_itr_t(); // set null
    }
  }
  v = ee->oppo()->vert();
  if (v->edge() == ee->oppo()) {
    m()[v].edge() = ee->prev();
    if (v->edge() == ee->oppo())
    {
      m()[v].edge() = edge_itr_t();
    }
  }

  // adjust the next and prev edges that the adjacent edges
  m()[ee->next()].prev() = ee->oppo()->prev();
  m()[ee->oppo()->prev()].next() = ee->next();
  m()[ee->prev()].next() = ee->oppo()->next();
  m()[ee->oppo()->next()].prev() = ee->prev();

  m().del(ee->oppo());
  m().del(ee);
  
  return 0;
}

//! @brief del face
//! @NOTE only delete face
template <typename TMPL>
int del(mesh_tmpl<TMPL> &m, const typename TMPL::faces_t::const_iterator& f)
{
  assert(f);
  using namespace std;
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  vector<edge_itr_t> fe;
  face_adj_edges<TMPL>(f, fe);

  const typename TMPL::faces_t::const_iterator ff = f;
  // set null face point
  for (auto eit = fe.begin(); eit != fe.end(); ++eit) {
    m()[(*eit)].face() = face_itr_t(); // null
    assert(!(*eit)->face());
  }

  assert(ff);
  m().del(ff);

  return 0;
}

//! @brief set opposite and boundary edge
//! used after adding all faces, the result is a halfedge mesh
template <typename TMPL>
int set_opposite_and_boundary_edge(mesh_tmpl<TMPL> &mt)
{
  using namespace std;
  typedef TMPL mesh_t;
  typedef const typename TMPL::vert_t* vert_ptr_t;
  typedef const typename TMPL::edge_t edge_t;
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  
  map<pair<vert_ptr_t, vert_ptr_t>, pair<edge_itr_t,edge_itr_t> > edges_map;

  mesh_t &m = mt();
  //edge_itr_t cei;
  edge_itr_t cei;
  for (cei=m.edges().begin(); cei!=m.edges().end(); ++cei) {
    pair<vert_ptr_t, vert_ptr_t> vptr;
    pair<edge_itr_t, edge_itr_t> eptr;

    bool is_swap = false;
    vptr.first = &(*(cei->prev()->vert()));
    vptr.second = &(*(cei->vert()));
    if (vptr.first == vptr.second) return 1; // a edge has two same vertices.
    if (vptr.first > vptr.second) {
      swap(vptr.first, vptr.second); // swap the addr
      is_swap = true;
    }
    auto em_it = edges_map.find(vptr);
    if (em_it == edges_map.end()) {  // the edge is the first time added to map
      if (is_swap) eptr.second = cei;
      else eptr.first = cei;
      edges_map.insert(make_pair(vptr,eptr));
    } else {
      if (!em_it->second.second) {
        if (is_swap) em_it->second.second = cei;
        else return 2; // the mesh isn't manifold, or the vertices sort of the corresponding face isn't right.
      } else if (!em_it->second.first) {
        if (!is_swap) em_it->second.first = cei;
        else return 3; // the mesh isn't manifold, or the vertices sort of the corresponding face isn't right.
      } else { return 4; }// the edge has more than two adjacent faces.
    }
  }

  map<vert_ptr_t, vector<edge_itr_t> > bound_vert;

  // add opposite relation
  for (auto em_it = edges_map.begin(); em_it != edges_map.end(); ++em_it) {
    pair<edge_itr_t, edge_itr_t>& eptr = em_it->second;
    if (eptr.first && eptr.second) {
      m[eptr.first].oppo() = eptr.second;
      m[eptr.second].oppo() = eptr.first;
    } else {
      edge_itr_t be = (!eptr.second)? eptr.first:eptr.second;
      assert(be);
      edge_itr_t ne = m.add(edge_t());   // add the new edge
      m[ne].vert() = be->prev()->vert(); // vert that new edge point to
      m[ne].oppo() = be;                 // add oppo of new edge
      m[be].oppo() = ne;                 // add oppo of current edge
      vert_ptr_t vptr = &(*(ne->vert()));
      auto bv_it = bound_vert.find(vptr);
      if (bv_it == bound_vert.end()) {
        std::vector<edge_itr_t> ve_itr;
        ve_itr.push_back(ne);
        bound_vert.insert(make_pair(vptr, ve_itr));
      } else {
        bv_it->second.push_back(ne);
      }
    }
  }

  // add next and prev relation of boundary edges.
  for (auto bv_it = bound_vert.begin(); bv_it != bound_vert.end(); ++bv_it) {
    vector<edge_itr_t> bound_edges;
    vector<edge_itr_t>& v_edges = bv_it->second; // in edges
    for (auto ve_it = v_edges.begin(); ve_it != v_edges.end(); ++ve_it) {
      edge_itr_t out_e = (*ve_it)->oppo();
      while (out_e->face()) { out_e = out_e->prev()->oppo(); }
      bound_edges.push_back(*ve_it); // in edge
      bound_edges.push_back(out_e);  // out edge
    }

    for (size_t i = 1; i+1 < bound_edges.size(); i+=2) {
      m[bound_edges[i]].prev() = bound_edges[i+1];
      m[bound_edges[i+1]].next() = bound_edges[i];
    }
    assert(bound_edges.size() > 1);
    m[bound_edges.front()].next() = bound_edges.back();
    m[bound_edges.back()].prev() = bound_edges.front();
  }

  // adjust the edge that verts point to
  for (auto vi = m().verts().begin(); vi != m.verts().end(); ++vi) {
    adjust_vert_edge(m, vi);
  }
  
  return 0;
}

//! @brief copy a half-edge mesh to another empty mesh
template <typename TMPL_A, typename TMPL_B>
int copy(const mesh_tmpl<TMPL_A>& ma, mesh_tmpl<TMPL_B>& mb)
{
  assert(mb().edges().size() == 0 && mb().faces().size() == 0
         && mb().verts().size() == 0);
  typedef typename TMPL_A::edge_t a_edge_t;
  typedef typename TMPL_B::edge_t b_edge_t;
  typedef typename TMPL_A::face_t a_face_t;
  typedef typename TMPL_B::face_t b_face_t;
  typedef typename TMPL_A::vert_t a_vert_t;
  typedef typename TMPL_B::vert_t b_vert_t;
  typedef typename TMPL_B::edges_t::const_iterator b_edge_itr_t;
  typedef typename TMPL_B::faces_t::const_iterator b_face_itr_t;
  typedef typename TMPL_B::verts_t::const_iterator b_vert_itr_t;
  
  std::map<const a_vert_t*, b_vert_itr_t> vert_map;
  std::map<const a_edge_t*, b_edge_itr_t> edge_map;
  std::map<const a_face_t*, b_face_itr_t> face_map;

  for (auto a_vi = ma().verts().begin(); a_vi != ma().verts().end(); ++a_vi) {
    b_vert_itr_t bv_itr = mb().add(b_vert_t());
    mb()[bv_itr].assign(*a_vi);
    vert_map.insert(std::make_pair(&(*a_vi), bv_itr));
  }

  for (auto a_ei = ma().edges().begin(); a_ei != ma().edges().end(); ++a_ei) {
    b_edge_itr_t be_itr = mb().add(b_edge_t());
    mb()[be_itr].assign(*a_ei);
    edge_map.insert(std::make_pair(&(*a_ei), be_itr));
  }

  for (auto a_fi = ma().faces().begin(); a_fi != ma().faces().end(); ++a_fi) {
    b_face_itr_t bf_itr = mb().add(b_face_t());
    mb()[bf_itr].assign(*a_fi);
    face_map.insert(std::make_pair(&(*a_fi), bf_itr));
  }

  for (auto a_vi = ma().verts().begin(); a_vi != ma().verts().end(); ++a_vi) {
    auto b_vi = vert_map[&(*a_vi)];
    mb()[b_vi].edge() = edge_map[&(*(a_vi->edge()))];
  }

  for (auto a_ei = ma().edges().begin(); a_ei != ma().edges().end(); ++a_ei) {
    auto b_ei = edge_map[&(*a_ei)];
    if (a_ei->face()) mb()[b_ei].face() = face_map[&(*(a_ei->face()))];
    mb()[b_ei].vert() = vert_map[&(*(a_ei->vert()))];
    mb()[b_ei].next() = edge_map[&(*(a_ei->next()))];
    mb()[b_ei].prev() = edge_map[&(*(a_ei->prev()))];
    mb()[b_ei].oppo() = edge_map[&(*(a_ei->oppo()))];
    assert(b_ei->vert());
    assert(b_ei->next());
    assert(b_ei->prev());
    assert(b_ei->oppo());
  }

  for (auto a_fi = ma().faces().begin(); a_fi != ma().faces().end(); ++a_fi) {
    auto b_fi = face_map[&(*a_fi)];
    mb()[b_fi].edge() = edge_map[&(*(a_fi->edge()))];
  }
  
  return 0;
}

//! @brief get the number of edges of face f
template<typename TMPL>
size_t valence(const typename TMPL::faces_t::const_iterator &f)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = f->edge();
  size_t r = 0;
  assert(e);
  do {
    ++r; e = e->next();
  } while (e != f->edge());
  return r;
}

//! @brief get the number of in edges of vert v
template<typename TMPL>
size_t valence(const typename TMPL::verts_t::const_iterator &v)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = v->edge();
  size_t r = 0;
  assert(e);
  do {
    ++r; e = e->next()->oppo();
  } while (e != v->edge());
  return r;
}

//! @brief get the adjacent edges of face f
template <typename TMPL, typename CON>
void face_adj_edges(const typename TMPL::faces_t::const_iterator &f, CON& fe)
{
  // assume the f is a valid iterator
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = f->edge();
  assert(e);
  do {
    fe.push_back(e);
    e = e->next();
  } while (e != fe.front());
}

//! @brief get the adjacent verts of face f
template <typename TMPL, typename CON>
void face_adj_verts(const typename TMPL::faces_t::const_iterator &f, CON& fv)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = f->edge();
  do {
    fv.push_back(e->vert());
    e = e->next();
  } while (e != f->edge());
}

//! @brief get the adjacent faces of face f
template <typename TMPL, typename CON>
void face_adj_faces(const typename TMPL::faces_t::const_iterator &f, CON& ff)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  edge_itr_t e = f->edge();
  do {
    face_itr_t f = e->oppo()->face();
    if (f) ff.push_back(f);
    e = e->next();
  } while (e != f->edge());
}

//! @brief get the adjacent faces of edge e
//! there are two faces adjacent an edge except the boundary edge
template <typename TMPL, typename CON>
void edge_adj_faces(const typename TMPL::edges_t::const_iterator &e, CON& ef)
{
  if (e->face()) ef.push_back(e->face());
  if (e->oppo()->face()) ef.push_back(e->oppo()->face());
}

//! @brief get the adjacent faces of vert v
template <typename TMPL, typename CON>
void vert_adj_faces(const typename TMPL::verts_t::const_iterator &v, CON& vf)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  edge_itr_t e = v->edge();
  do {
    face_itr_t f = e->face();
    if (f) vf.push_back(f);
    e = e->next()->oppo();
  } while (e != v->edge());
}

//! @brief get the adjacent verts of vert v
template <typename TMPL, typename CON>
void vert_adj_verts(const typename TMPL::verts_t::const_iterator &v, CON& vv)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;
  edge_itr_t e = v->edge();
  do {
    vert_itr_t v = e->vert();
    assert(v);
    vv.push_back(v);
    e = e->next()->oppo();
  } while (e != v->edge());
}

//! @brief get the adjacent out edges of vert v
template <typename TMPL, typename CON>
void vert_adj_out_edges(const typename TMPL::verts_t::const_iterator &v, CON& ve)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = v->edge()->oppo();
  do {
    ve.push_back(e);
    e = e->prev()->oppo();
  } while (e != v->edge()->oppo());
}

//! @brief get the half edge of two vertices.
template <typename TMPL>
typename TMPL::edges_t::const_iterator
get_edge(const typename TMPL::verts_t::const_iterator & v1,
         const typename TMPL::verts_t::const_iterator & v2)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = v2->edge();
  do
  {
    if (!e) break;
    if (e->oppo()->vert() == v1) return e;
    e = e->next()->oppo();
  } while (e != v2->edge());
  return edge_itr_t();
}


template <typename TMPL, typename CVI>
void
adjust_vert_edge(mesh_tmpl<TMPL> &m, const CVI & vert)
{
  //return for condition already holds
  if (!vert->edge() || !vert->edge()->face())
  {
    return;
  }
  
  typedef typename TMPL::edges_t::const_iterator CEI;
  CEI ei = vert->edge();

  while (ei)
  {
    if (!(ei->face()) )
    {
      m()[vert].edge() = ei;
      break;
    }
    ei = ei->next()->oppo();
    
    if (ei == vert->edge())
    {
      break;
    }
  }
}

// @brief all edges in sec are border edge, and the region between
// sec[i*2] and sec[i*2+1] are connected by faces.
template <typename CVI, typename Con>
void sectors(const CVI &vi, Con &sec)
{
  auto ei = vi->edge();
  do
  {
    if(!ei->face())
    {
      //add in and out edges of a sector sequencely
      sec.push_back(ei);
      sec.push_back(ei->next());
    }
    ei = ei->next()->oppo();
  }
  while (ei != vi->edge());
}


//! @brief add a new sector from in to out into sec., in and out are
//! outer edge of a sector, but may be not border edges.
template <typename TMPL, typename CEI, typename Con>
int add_face_into_sectors(mesh_tmpl<TMPL> &m, Con &sec, const CEI &in, const CEI &out)
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
add_face2(mesh_tmpl<TMPL> &m, const RND_ITR &vert_loop_beg, std::size_t n)
{
  //add vertices
  std::vector<bool> is_v_exist(n, true);//is the vertex already exists
  for (size_t i=0; i<n; ++i)
  {
    is_v_exist[i] =!!(vert_loop_beg[i]->edge());
  }
  
  assert(n > 2); // TODO: may be also good for n=1,2.


  typedef typename TMPL::faces_t::const_iterator CFI;
  CFI fi = m().add(typename TMPL::face_t());
  typedef typename TMPL::edges_t::const_iterator CEI;
  typedef std::vector<CEI> e_con;
  e_con edges(n); // edges[i]: v_i -> v_{i+1}

  // add edges, set vert and oppo
  for(size_t i = 0; i < n; ++i) {
    // find edge v_i to v_{i+1}
    CEI ei = get_edge<TMPL>(vert_loop_beg[i], vert_loop_beg[(i+1)%n]);
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
      std::vector<CEI> bd;
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
    if(!vert_loop_beg[i]->edge())
    {
      //opposite edge is more likely a boundary edge,
      //for effciency, not use edges[i-1]
      m()[vert_loop_beg[i]].edge() = edges[i]->oppo();
    }
    m()[edges[i]].face() = fi; // the edge adjacent to this face
    assert(vert_loop_beg[i]->edge()->vert() == vert_loop_beg[i]);
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


//! @brief is an given halfedge can be flipped
//! @param is_valid(i) = true
//! @return 0 for done, 1 for incorrect edge denoting, 2 for non-triangle
//! 3 for special cases
template <typename TMPL>
int try_edge_flip(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;
  //get faces adjacent to edge
  edge_itr_t e_op = e->oppo();
  
  //whether there is two face adjacent to the given edge
  if (!e->face() || !e_op->face())
  {
    return 1;
  }
  
  //whether faces adjacent to the edge are triangles
  if (valence<TMPL>(e->face()) != 3
      ||valence<TMPL>(e_op->face()) != 3)
  {
    return 2;
  }

  //special cases can not flip
  //there is an edge in the position which we want to add new edge
  auto sv = e->next()->vert();
  auto ev = e->oppo()->next()->vert();
  if (!!get_edge<TMPL>(sv, ev))
  {
    return 3;
  }
  
  return 0;
}


//! @brief implement a flip operation on a given edge by rotate the
//! existing faces
//! @param m for mesh, e for the edge to be flipped
//! @warning {faces are preserved with some of their member being changed
//! so this may cause some problems when deal with fetures of these faces}
template <typename TMPL>
void edge_flip_by_rotate(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;

  //get opposite half_edge
  edge_itr_t e_op = e->oppo();
  
  //flip
  //modify vertices
  m()[e].vert() = e->next()->vert();
  m()[e_op].vert() = e_op->next()->vert();
  if (e->vert()->edge() == e)
  {
    //modify vertex whose edge is e_op
    m()[e->vert()].edge() = e->next()->oppo();
  }
  if (e_op->vert()->edge() == e_op)
  {
    //modify vertex whose edge is e
    m()[e_op->vert()].edge() = e_op->next()->oppo();
  }

  //modify faces
  m()[e->face()].edge() = e;
  m()[e_op->face()].edge() = e_op;

  //modify edges
  m()[e->next()].prev() = e_op->prev();
  m()[e->next()].next() = e_op;
  m()[e->next()].face() = e_op->face();
  
  m()[e->prev()].next() = e_op->next();
  m()[e->prev()].prev() = e;

  m()[e_op->next()].prev() = e->prev();
  m()[e_op->next()].next() = e;
  m()[e_op->next()].face() = e->face();
  
  m()[e_op->prev()].next() = e->next();
  m()[e_op->prev()].prev() = e_op;

  //modify e and e_op, use tmp space to save pointer
  edge_itr_t tmp = e->next();
  m()[e].next() = e->prev();
  m()[e].prev() = e_op->next();

  m()[e_op].next() = e_op->prev();
  m()[e_op].prev() = tmp;
}

//! @param try_edge_flip(m, i) == 0;
template <typename TMPL>
void edge_flip_by_del_add(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;

  std::vector<vert_itr_t> new_face, new_face2;
  new_face.push_back(e->vert());
  new_face.push_back(e->next()->vert());
  new_face.push_back(e->oppo()->next()->vert());
  new_face2.push_back(e->oppo()->vert());
  new_face2.push_back(e->oppo()->next()->vert());
  new_face2.push_back(e->next()->vert());
  del(m, e);
  add_face2(m, new_face.begin(), 3);
  add_face2(m, new_face2.begin(), 3);
}

//! @brief test whether an edge can be collapse
template <typename TMPL>
int try_collapse(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;

  std::vector<edge_itr_t> e_svert;
  std::vector<edge_itr_t> e_evert;
  
  vert_adj_out_edges<TMPL>(e->oppo()->vert(), e_svert);
  vert_adj_out_edges<TMPL>(e->vert(), e_evert);

  auto sitr = e_svert.begin();
  auto eitr = e_evert.begin();

  for (; sitr!=e_svert.end(); ++sitr)
  {
    if (*sitr == e->prev()->oppo() || *sitr == e->oppo()->next())
    {
      continue;
    }
    for (; eitr!=e_evert.end(); ++eitr)
    {
      //if there are more than two faces share the collapsing edge
      if (*eitr == e->next() || *eitr == e->oppo()->prev()->oppo())
      {
        continue;
      }
      if ((*eitr)->vert() == (*sitr)->vert())
      {
        return 2;
      }
    }
  }
  
  return 0;
}

//! @brief implement a collapse operation on a given halfedge
template <typename TMPL>
const typename TMPL::verts_t::const_iterator &
collapse_edge(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  const typename TMPL::verts_t::const_iterator& rtn = e->vert();
  
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  
  //save e->oppo() and start vertex before it maybe invisable later
  auto e_op = e->oppo();
  auto e_vert = e->vert(), s_vert = e_op->vert();
  
  //special cases
  //the last edges
  if (e->next() == e->prev())
  {
    m()[e->vert()].edge() = edge_itr_t();//set to null
    //delete the collapsed edge
    m().del(e);
    m().del(e_op);
    //delete the start vertex
    m().del(s_vert);
    return rtn;
  }

  //NOTE:1. e->vert()->edge() should not be used in the remaining
  //part of this function for it had been changed
  //2. can not move this part to the end because if topology information
  //changes, we may fail to find vert->edge();
  //3. a basic thought, not change topology info unless it's necessary
  
  //modify e->vert->edge() for two condition
  //#1 e->oppo()->prev() will be deleted if valence of that face is 3
  if (e->oppo()->face()
      && valence<TMPL>(e->oppo()->face()) == 3
      && e->vert()->edge() == e->oppo()->prev())
  {
    m()[e->vert()].edge() = e->oppo()->next()->oppo();
  }

  //#2 
  if (e->vert()->edge() == e)
  {
    if (e->face())
    {
      m()[e->vert()].edge() = e->next()->oppo();
    }
    else
    {
      if (e->oppo()->face())
      {
        m()[e->vert()].edge() = e->oppo()->next()->oppo();
      }
      else
      {
        //in this else scope, both edge have no faces
        if (e->next() != e->oppo())
        {
          m()[e->vert()].edge() = e->next()->oppo();
        }
        else
        {
          if (e->prev() != e->oppo())
          {
            m()[e->vert()].edge() = e->prev();
          }
          else
          {
            m()[e->vert()].edge() = edge_itr_t();
          }
        }
      }
    }
  }

  
  //normal procedure
  //repoint all the edges to the end vertex of e
  std::vector<edge_itr_t> edges_adj_svert;
  vert_adj_out_edges<TMPL>(s_vert, edges_adj_svert);
  
  auto eitr = edges_adj_svert.begin();
  for (;eitr!=edges_adj_svert.end(); ++eitr)
  {
    if (*eitr == e) continue;

    //delete loop
    if ((*eitr)->vert() == e_vert())
    {
      if ((*eitr)->vert()->edge() == (*eitr)->oppo())
      {
        m()[(*eitr)->vert()].edge() = (*eitr)->next();
      }
      if ((*eitr)->oppo()->vert()->edge() == (*eitr))
      {
        m()[(*eitr)->oppo()->vert()].edge() = (*eitr)->oppo()->next();
      }
      m()[(*eitr)->next()].prev() = (*eitr)->prev();
      m()[(*eitr)->prev()].next() = (*eitr)->next();
      m()[(*eitr)->oppo()->next()].prev() = (*eitr)->oppo()->prev();
      m()[(*eitr)->oppo()->prev()].next() = (*eitr)->oppo()->next();
      m().del((*eitr)->oppo());
      m().del(*eitr);
    }
    //repoint
    m()[(*eitr)->oppo()].vert() = e_vert;
  }

  
  //modify edges
  if (e->face())
  {
    if (valence<TMPL>(e->face()) == 3)
    {
      
      m()[e->prev()->oppo()].oppo() = e->next()->oppo();
      m()[e->next()->oppo()].oppo() = e->prev()->oppo();
      if (e->next()->vert()->edge() == e->next())
      {
        m()[e->next()->vert()].edge() = e->prev()->oppo();
      }
      m().del(e->next());
      m().del(e->prev());
      m().del(e->face());
    }
    else
    {
      //faces with more than 3 edges
      if (e->face()->edge() == e)
      {
        m()[e->face()].edge() = e->next();
      }
      m()[e->prev()].next() = e->next();
      m()[e->next()].prev() = e->prev();
    }
    //if edges with no face need to be delete, add it here
  }
  else
  {
    //NOTE: e->prev() or e->next() may equal to e_op
    if (e->next() == e_op)
      m()[e->prev()].next() = e_op->next();
    else
      m()[e->prev()].next() = e->next();

    if (e->prev() == e_op)
      m()[e->next()].prev() = e_op->prev();
    else
      m()[e->next()].prev() = e->prev();
  }

  
  if (e_op->face())
  {
    if (valence<TMPL>(e_op->face()) == 3)
    {
      m()[e_op->next()->oppo()].oppo() = e_op->prev()->oppo();
      m()[e_op->prev()->oppo()].oppo() = e_op->next()->oppo();

      //modify edge of vertex before delete edge
      if (e_op->next()->vert()->edge() == e_op->next())
      {
        m()[e_op->next()->vert()].edge() = e_op->prev()->oppo();
      }
      m().del(e_op->next());
      m().del(e_op->prev());
      m().del(e_op->face());
    }
    else
    {
      if (e_op->face()->edge() == e_op)
      {
        m()[e_op->face()].edge() = e_op->next();
      }
      
      m()[e_op->prev()].next() = e_op->next();
      m()[e_op->next()].prev() = e_op->prev();
    }
  }
  else
  {
    //NOTE: e_op->prev() or e_op->next() may equal to e
    if (e_op->next() == e)
      m()[e_op->prev()].next() = e->next();
    else
      m()[e_op->prev()].next() = e_op->next();

    if (e_op->prev() == e)
      m()[e_op->next()].prev() = e->prev();
    else
      m()[e_op->next()].prev() = e_op->prev();
  }

  //delete the collapsed edge
  m().del(e);
  m().del(e_op);
  //delete the start vertex
  m().del(s_vert);

  adjust_vert_edge(m, rtn);
  return rtn;
}

template <typename TMPL>
int is_collapse_ok(mesh_tmpl<TMPL> &m, const typename TMPL::verts_t::const_iterator &v)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  
  std::vector<edge_itr_t> ve;
  vert_adj_out_edges<TMPL>(v, ve);

  //merge two edges which link same vertices
  //compare from the last edge
  auto last_e = ve.end();
  --last_e;
  auto itr = ve.begin();

  while (itr != ve.end())
  {
    if (!(*itr)->face() && !(*itr)->oppo()->face())
    {
      collapse_edge(m, *itr);
      ++itr;
      continue;
    }

    //NOTE this part is prepare for non-manifold cases
    if ((*itr)->vert() == (*last_e)->vert())
    {
      //merge two adjcent edges together
      if ((*itr)->face() || (*last_e)->oppo()->face())
      {
        //case1
        if (v->edge() == (*last_e))
        {
          m()[v].edge() = *itr;
        }
        if ((*last_e)->vert()->edge() == (*itr)->oppo())
        {
          m()[(*last_e)->vert()].edge() = (*last_e)->oppo();
        }

        m()[(*last_e)->oppo()].oppo() = *itr;
        m()[(*itr)].oppo() = (*last_e)->oppo();
        m().del(*last_e);
        m().del((*itr)->oppo());
      }
      else
      {
        //case2
        if (v->edge() == *itr)
        {
          m()[v].edge() = *last_e;
        }
        if ((*itr)->vert()->edge() == (*last_e)->oppo())
        {
          m()[(*itr)->vert()].edge() = (*itr)->oppo();
        }
        
        m()[(*last_e)].oppo() = (*itr)->oppo();
        m()[(*itr)->oppo()].oppo() = *last_e;
        m().del((*last_e)->oppo());
        m().del(*itr);
      }
    }
  
    last_e = itr;
    ++itr;
  }
  return 0;
}

//! @brief split 1 edge into 2 edges, 2 faces into 4 faces and fix 
//! the adj of relational verts, edges, faces.
//! the adjacent faces of the edge should be triangles.
template <typename TMPL>
const typename TMPL::verts_t::const_iterator
split_edge(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;  
  if (!e) return vert_itr_t();

  const size_t is_triangle = 3;

  if ((e->face()) && valence<TMPL>(e->face())!=is_triangle)
    return vert_itr_t();
  if ((e->oppo()->face()) && valence<TMPL>(e->oppo()->face())!=is_triangle)
    return vert_itr_t();

  vert_itr_t new_vert = m().add(typename TMPL::vert_t());
  m()[new_vert].edge() = e->oppo();

  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  
  edge_itr_t mid_e1 = m().add(typename TMPL::edge_t());
  edge_itr_t mid_e2 = m().add(typename TMPL::edge_t());
  
  m()[mid_e1].vert() = new_vert; m()[mid_e1].oppo() = mid_e2;
  m()[mid_e1].next() = e; m()[mid_e1].prev() = e->prev();
  m()[mid_e1].face() = face_itr_t();

  m()[e->prev()].next() = mid_e1;
  
  m()[mid_e2].vert() = e->oppo()->vert();
  m()[mid_e2].oppo() = mid_e1; m()[mid_e2].next() = e->oppo()->next();  
  m()[mid_e2].prev() = e->oppo(); m()[mid_e2].face() = face_itr_t();

  m()[e->oppo()->next()].prev() = mid_e2;

  if (e->oppo()->vert()->edge() == e->oppo())
    m()[e->oppo()->vert()].edge() = mid_e2;

  m()[e].prev() = mid_e1;
  m()[e->oppo()].next() = mid_e2;
  m()[e->oppo()].vert() = new_vert;

  // split left face
  if (e->face()) {
    m()[e->face()].edge() = e;
    
    face_itr_t left_f = m().add(typename TMPL::face_t());
    m()[left_f].edge() = mid_e1;
    
    m()[mid_e1].face() = left_f; m()[mid_e1->prev()].face() = left_f;
    
    edge_itr_t left_e1 = m().add(typename TMPL::edge_t());
    edge_itr_t left_e2 = m().add(typename TMPL::edge_t());
    
    m()[left_e1].vert() = new_vert;
    m()[left_e1].next() = e; m()[left_e1].prev() = e->next();
    m()[left_e1].oppo() = left_e2; m()[left_e1].face() = e->face();
    
    m()[left_e2].vert() = e->next()->vert();
    m()[left_e2].next() = mid_e1->prev(); m()[left_e2].prev() = mid_e1;
    m()[left_e2].oppo() = left_e1; m()[left_e2].face() = left_f;
    
    m()[e].prev() = left_e1; m()[e->next()].next() = left_e1;
    m()[mid_e1].next() = left_e2; m()[mid_e1->prev()].prev() = left_e2;
  }

  // split right face
  if (e->oppo()->face()) {
    m()[e->oppo()->face()].edge() = e->oppo();
    
    face_itr_t right_f = m().add(typename TMPL::face_t());
    m()[right_f].edge() = mid_e2;

    m()[mid_e2].face() = right_f; m()[mid_e2->next()].face() = right_f;

    edge_itr_t right_e1 = m().add(typename TMPL::edge_t());
    edge_itr_t right_e2 = m().add(typename TMPL::edge_t());

    m()[right_e1].vert() = mid_e2->next()->vert();
    m()[right_e1].next() = e->oppo()->prev();
    m()[right_e1].prev() = e->oppo();
    m()[right_e1].oppo() = right_e2;
    m()[right_e1].face() = e->oppo()->face();
    m()[right_e2].vert() = new_vert;
    m()[right_e2].next() = mid_e2;
    m()[right_e2].prev() = mid_e2->next();
    m()[right_e2].oppo() = right_e1;
    m()[right_e2].face() = right_f;

    m()[e->oppo()->prev()].prev() = right_e1;
    m()[mid_e2->next()].next() = right_e2;
    m()[e->oppo()].next() = right_e1;
    m()[mid_e2].prev() = right_e2;
  }

  if (!e->face()) m()[new_vert].edge() = mid_e1;

  return new_vert;
}
