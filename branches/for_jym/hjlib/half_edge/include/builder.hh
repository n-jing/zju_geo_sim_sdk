#ifdef ZJUCAD_MATRIX_VERSION

template <typename ME>
mat_builder<ME>
build(const zjucad::matrix::matrix_expression<ME> &cells) {
  return mat_builder<ME>(cells());
}

template <typename TMPL, typename ME>
void
assign(mesh_tmpl<TMPL> &mt, const mat_builder<ME> &b)
{
  //assign_from_mat_with_addface(mt, b.cells_);
  assign_from_mat(mt, b.cells_);
}

//! @assume triangular manifold
template <typename TMPL, typename ME>
void
assign_from_mat(mesh_tmpl<TMPL> &mt, const zjucad::matrix::matrix_expression<ME> &cells)
{
  assert(cells().size(1) == 3);
  using namespace zjucad::matrix;
  using namespace std;

  typedef TMPL mesh_t;
  typedef typename mesh_t::verts_t::const_iterator CVI;
  typedef typename mesh_t::edges_t::const_iterator CEI;
  mesh_t &m = mt();

  // insert vert
  const size_t vert_num = zjucad::matrix::max(cells())+1;
  matrix<CVI> vis(vert_num);
  for(size_t vi = 0; vi < vert_num; ++vi)
    vis[vi] = m.add(typename mesh_t::vert_t());

  // insert (inner) edges and face without set oppo
  vector<map<size_t, CEI> > adj_map(vert_num); // adj_map[from][to] = he
  for(size_t fi = 0; fi < cells().size(2); ++fi) {
    typename TMPL::faces_t::const_iterator f = m().add(typename TMPL::face_t());
    typedef std::vector<typename TMPL::edges_t::const_iterator> e_con;
    e_con edges;
    edges.reserve(cells().size(1));

    // insert
    for(size_t ei = 0; ei < cells().size(1); ++ei)
      edges.push_back(m().add(typename TMPL::edge_t()));
    // set edge and vert, init map
    for(size_t ei = 0; ei < cells().size(1); ++ei) {
      if(!vis[cells()(ei, fi)]->edge())
        m[vis[cells()(ei, fi)]].edge() = edges[ei];

      m[edges[ei]].vert() = vis[cells()(ei, fi)];
      m[edges[ei]].face() = f;
      m[edges[ei]].next() = edges[(ei+1)%edges.size()];
      m[edges[ei]].prev() = edges[(edges.size()+ei-1)%edges.size()];

      const size_t from = cells()(ei, fi), to = cells()((ei+1)%cells().size(1), fi);
      auto loc = adj_map[from].find(to);
      if(loc != adj_map[from].end()) {
        char error[256];
        sprintf(error, "in assign_from_mat: duplicated edge from %ld to %ld.", from, to);
        throw std::runtime_error(error);
      }
      else
        adj_map[from].insert(loc, make_pair(to, edges[(ei+1)%cells().size(1)]));
    }
    // set face
    m[f].edge() = edges[0];
  }

  // add boundary edges and set oppo
  for(size_t vi = 0; vi < adj_map.size(); ++vi) {
    for(typename map<size_t, CEI>::const_iterator vjj = adj_map[vi].begin();
        vjj != adj_map[vi].end(); ++vjj) {
      const size_t vj = vjj->first;
      assert(vi != vj);
      CEI &oppo = adj_map[vj][vi];
      if(!oppo) { // pair the boundary
        oppo = m().add(typename TMPL::edge_t());
        m[oppo].vert() = vjj->second->prev()->vert();
      }
      m[vjj->second].oppo() = oppo;
      m[oppo].oppo() = vjj->second;
    }
  }

  // link boundary edges
  for(CEI ei = m.edges().begin(); ei != m.edges().end(); ++ei) {
    if(ei->face()) {
      assert(!!ei->next() && !!ei->prev());
      continue;
    }
    for(CEI ej = ei->oppo();;) {
      ej = ej->prev()->oppo(); // ccw of out edges
      if(!ej->face()) { // border
        m[ei].next() = ej;
        m[ej].prev() = ei;
        break;
      }
    }
  }
}

//! @assign from mat use add_face function
//! @NOTE use add_face function will make each halfedge and its
//! @opposite one adjacent in a vector thus improve cache hit proportion
template <typename TMPL, typename ME>
void
assign_from_mat_with_addface(mesh_tmpl<TMPL> &mt, const zjucad::matrix::matrix_expression<ME> &cells)
{
  assert(cells().size(1) == 3);
  using namespace zjucad::matrix;
  using namespace std;
  typedef TMPL mesh_t;
  typedef typename mesh_t::verts_t::const_iterator CVI;
  typedef typename mesh_t::edges_t::const_iterator CEI;
  mesh_t &m = mt();

  // insert vert
  const size_t vert_num = zjucad::matrix::max(cells())+1;
  vector<CVI> vis;
  vis.reserve(vert_num);
  for(size_t vi = 0; vi < vert_num; ++vi)
    vis[vi] = m.add(typename mesh_t::vert_t());

  for (size_t fi=0; fi<cells().size(2); ++fi)
  {
    vector<CVI> vert_loop;
    vert_loop.clear();
    vert_loop.reserve(3);
    for (size_t i=0; i<cells().size(1); ++i)
    {
      vert_loop.push_back(vis[cells()(i, fi)]);
    }

    add_face2(mt, vert_loop.begin(), cells().size(1));
  }
}

//! @assume cells code the face list of a triangular manifold
template <typename ME>
mat_builder_with_id<ME>
build_with_id(const zjucad::matrix::matrix_expression<ME> &cells)
{
  return mat_builder_with_id<ME>(cells());
}

template <typename TMPL, typename ME>
void
assign(mesh_tmpl<TMPL> &mt, const mat_builder_with_id<ME> &b)
{
  assign_from_mat(mt, b.cells_);

  typedef TMPL mesh_t;

  typedef typename mesh_t::verts_t::const_iterator CVI;
  typedef typename mesh_t::faces_t::const_iterator CFI;

  size_t id = 0;
  for(auto vi = mt().verts().begin(); vi != mt().verts().end(); ++vi, ++id)
    mt()[vi].id() = id;

  id = 0;
  for(auto fi = mt().faces().begin(); fi != mt().faces().end(); ++fi, ++id)
    mt()[fi].id() = id;
}

#endif
