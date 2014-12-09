

template <typename TMPL>
bool read_mesh(TMPL & m, const char* _filename)
{

  typedef zjucad::matrix::matrix<size_t> matrixst;
  zjucad::matrix::matrix<size_t> tris;
  zjucad::matrix::matrix<double> nods;
  if (jtf::mesh::load_obj(_filename, tris, nods))
  {
    std::cout << "load obj failed!" << std::endl;
    exit(0);
  }
  m = build(tris);
  size_t idx = 0;
  for(auto i = m.verts().begin(); i != m.verts().end(); ++i, ++idx)
    assign(m[i].coord_, nods(colon(), idx));
  return true;
}

template <typename TMPL>
bool write_mesh(TMPL & m, const char* _filename)
{
  typedef TMPL mesh_t;
  typedef typename mesh_t::verts_t::const_iterator CVI;
  typedef typename mesh_t::faces_t::const_iterator CFI;
  typedef typename mesh_t::edges_t::const_iterator CEI;

  using namespace hj::half_edge;
  using namespace zjucad::matrix;
  
  typedef zjucad::matrix::matrix<size_t> matrixst;
  
  matrix<size_t> tris(3, m.faces().size());
  matrix<double> nods(3, m.verts().size());

  CVI vi=m.verts().begin();
  for (size_t i=0; vi!=m.verts().end(); ++vi)
  {
    for (size_t j=0; j<3; ++j)
      nods(j, i) = vi->coord_(j, 0);
    //assign(nods(colon(),i), vi->coord_);
    m()[vi].id_ = i++;
  }

  CFI fi = m.faces().begin();
  for (size_t j=0; fi!=m.faces().end(); ++fi, ++j)
  {
    CEI ei = fi->edge();
    for (size_t i=0; i<3; ++i)
    {
      tris(i, j) = ei->vert()->id_;
      ei = ei->next();
    }
  }

  if(jtf::mesh::save_obj(_filename, tris, nods))
    return true;
  else
    return false;
}
