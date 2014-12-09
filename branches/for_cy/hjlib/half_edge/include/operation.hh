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
    m()[ei].face() = fi; // the edge adjacent to this face
  }
  vi = vert_loop_beg;
  for(std::size_t i = 0; i < edges.size(); ++i, ++vi) {
    if(!(*vi)->edge())
      m()[*vi].edge() = edges[(i+1)%edges.size()]; // ei points to i
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
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = v->edge();
  if (!e) return false; // the vert is isolated vertex
  do {
    if (e && is_boundary(e)) return true;
    e = e->next();
  } while (e != v->edge());
  
  return true;
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
  
  bool is_del = !(ee->oppo()->face());
  if (ee->face()) {
    if (del(m, ee->face())) return 1;
  }
  if (!is_del) { // can't use e->oppo()->face() directly, because the e may be deleted
    if (del(m, ee->oppo()->face())) return 1;
  }

  // adjust the edge that vert point to
  vert_itr_t v = ee->vert();
  if (v->edge() == ee->oppo()) {
    m()[v].edge() = ee->next();
    if (v->edge() == ee->oppo()) m()[v].edge() = edge_itr_t(); // set null
  }
  v = ee->oppo()->vert();
  if (v->edge() == ee) {
    m()[v].edge() = ee->oppo()->next();
    if (v->edge() == ee) m()[v].edge() = edge_itr_t();
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

//! @brief get the number of out edges of vert v
template<typename TMPL>
size_t valence(const typename TMPL::verts_t::const_iterator &v)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = v->edge();
  size_t r = 0;
  assert(e);
  do {
    ++r; e = e->prev()->oppo();
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
    e = e->prev()->oppo();
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
    e = e->prev()->oppo();
  } while (e != v->edge());
}

//! @brief get the adjacent out edges of vert v
template <typename TMPL, typename CON>
void vert_adj_out_edges(const typename TMPL::verts_t::const_iterator &v, CON& ve)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = v->edge();
  do {
    ve.push_back(e);
    e = e->prev()->oppo();
  } while (e != v->edge());
}

//! @brief get the half edge of two vertices.
template <typename TMPL>
typename TMPL::edges_t::const_iterator
get_edge(const typename TMPL::verts_t::const_iterator& v1,
         const typename TMPL::verts_t::const_iterator& v2)
{
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  edge_itr_t e = v1->edge();
  do {
    if (e->vert() == v2) return e;
    e = e->prev()->oppo();
  } while (e != v1->edge());
  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  return edge_itr_t();
}

//! TODO
template <typename TMPL>
int flip_edge(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  assert(1);
  return 0;
}

//! TODO
template <typename TMPL>
const typename TMPL::verts_t::const_iterator&
collapse_edge(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
  assert(0);
  return 0;
}


//! @brief split 1 edge into 2 edges, 2 faces into 4 faces and fix the adj of relational verts, edges, faces.
template <typename TMPL>
const typename TMPL::verts_t::const_iterator
split_edge(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e)
{
    assert(e);
    /**
     *                   e8|^
     *   |^            f2  ||   original face(of)         
     *   ||                V|e1                
     * e5||e1   ==> e7<---    <---e2                 
     *   ||         e6--->    --->e3                
     *   v|              e5|^             
     *             (of)    ||   f1        
     *                     V|e4
     */

  typedef typename TMPL::edges_t::const_iterator edge_itr_t;
  typedef typename TMPL::faces_t::const_iterator face_itr_t;
  typedef typename TMPL::verts_t::const_iterator vert_itr_t;

  // is triangle?
  const size_t my_Triangle = 3;
  assert(e->face());
  assert(e->oppo()->face());
  if ( (valence<TMPL>(e->face()) != my_Triangle) || valence<TMPL>(e->oppo()->face()) != my_Triangle )
      return vert_itr_t();
    
  vert_itr_t new_vert = m().add(typename TMPL::vert_t());
  m()[new_vert].edge() = e;

  edge_itr_t e1 = e;
  edge_itr_t e2 = m().add(typename TMPL::edge_t());
  edge_itr_t e3 = m().add(typename TMPL::edge_t());
  edge_itr_t e4 = m().add(typename TMPL::edge_t());

  edge_itr_t e5 = e->oppo();
  edge_itr_t e6 = m().add(typename TMPL::edge_t());
  edge_itr_t e7 = m().add(typename TMPL::edge_t());
  edge_itr_t e8 = m().add(typename TMPL::edge_t());

  face_itr_t f1 = m().add(typename TMPL::face_t());
  face_itr_t f2 = m().add(typename TMPL::face_t());


  // direct adj face
  // find the mid edge of a face
  edge_itr_t edge_mid = e1->next();
  vert_itr_t vert_mid = edge_mid->vert();
  // two faces
  m()[f1].edge() = e3;
  m()[e1->face()].edge() = e1;
  // four edges
  m()[e2].vert() = new_vert;  m()[e2].face() = e1->face();  m()[e2].prev() = edge_mid;   m()[e2].next() = e1;               m()[e2].oppo() = e3;
  m()[e3].vert() = vert_mid;  m()[e3].face() = f1;          m()[e3].prev() = e4;         m()[e3].next() = edge_mid->next(); m()[e3].oppo() = e2;
  m()[e4].vert() = new_vert;  m()[e4].face() = f1;          m()[e4].prev() = e1->prev(); m()[e4].next() = e3;               m()[e4].oppo() = e5;                              
  m()[e1].prev() = e2;                                            m()[e1].oppo() = e8;
  // other edges
  m()[e1->prev()].next() = e4;
  m()[edge_mid->next()].prev() = e3;
  m()[edge_mid].next() = e2;


  // oppo face
  edge_itr_t edge_oppo_mid = e5->next();
  vert_itr_t vert_oppo_mid = edge_oppo_mid->vert();
  // two faces
  m()[f2].edge() = e7;
  m()[e5->face()].edge() = e5;
  // four edges
  m()[e6].vert() = new_vert;      m()[e6].face() = e5->face(); m()[e6].prev() = edge_oppo_mid; m()[e6].next() = e5;                    m()[e6].oppo() = e7;
  m()[e7].vert() = vert_oppo_mid; m()[e7].face() = f2;         m()[e7].prev() = e8;            m()[e7].next() = edge_oppo_mid->next(); m()[e7].oppo() = e6;
  m()[e8].vert() = new_vert;      m()[e8].face() = f2;         m()[e8].prev() = e5->prev();    m()[e8].next() = e7;                    m()[e8].oppo() = e1;
  m()[e5].prev() = e6;                                                    m()[e5].oppo() = e4;
  // other edges
  m()[e5->prev()].next() = e8;
  m()[edge_oppo_mid->next()].prev() = e7;
  m()[edge_oppo_mid].next() = e6;


  return m()[e2].vert();
}

