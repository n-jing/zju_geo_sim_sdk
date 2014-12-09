#ifndef HJ_HALF_EDGE_OPERATION_H_
#define HJ_HALF_EDGE_OPERATION_H_

//! @file operation.h
//! @brief declaration of the basic operations on half_edge
//! @author Jin Huang, Xianzhong Fang
//! @date 201406-

#include <vector>
#include <cassert>
#include <map>
#include <algorithm>
#include <iostream>
#include "half_edge.h"
#include "container.h"

namespace hj { namespace half_edge {

//@{
//! @brief query the valid status of entity
template <typename TMPL>
bool
is_valid(const cont_tmpl<TMPL> &c, const typename TMPL::const_iterator &i);

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m, const typename TMPL::verts_t::const_iterator &vi);

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &ei);

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m, const typename TMPL::faces_t::const_iterator &fi);

template <typename TMPL, typename VC, typename EC, typename FC>
void
find_invalid(const mesh_tmpl<TMPL> &m, VC &vc, EC &ec, FC &fc);

template <typename TMPL>
bool
is_valid(const mesh_tmpl<TMPL> &m);
//@}

//@{
//! @brief judge if the element is a boundary element
//! boundary face: at least one adjacent edge's opposite edge is boundary.
//! boundary edge: the face that edge points to is null.
//! boundary vert: at least one adjacent edge is boundary.
template <typename TMPL>
    bool is_boundary(const typename TMPL::verts_t::const_iterator& v);
template <typename TMPL>
    bool is_boundary(const typename TMPL::edges_t::const_iterator& e);
template <typename TMPL>
    bool is_boundary(const typename TMPL::faces_t::const_iterator& f);
//@}

//@{
//! @brief judge if the input element is a isolated element
template <typename TMPL>
    bool is_isolated(const typename TMPL::verts_t::const_iterator& v);
template <typename TMPL>
    bool is_isolated(const typename TMPL::edges_t::const_iterator& e);
template <typename TMPL>
    bool is_isolated(const typename TMPL::faces_t::const_iterator& f);
//@}

//! @brief not automatically glue into the rest of the faces for
//! efficiency, leaving the mesh with 2-side holes.
template <typename TMPL, typename FWD_ITR>
typename TMPL::faces_t::const_iterator
add_face(mesh_tmpl<TMPL> &m, const FWD_ITR &vert_loop_beg, std::size_t n);

//! @brief NOTICE append to border_edges.
template <typename TMPL, typename C>
void
find_border(const cont_tmpl<TMPL> &edges, C &border_edges);

//! @brief del vert
template <typename TMPL>
    int del(mesh_tmpl<TMPL> &m, const typename TMPL::verts_t::const_iterator& v);

//! @brief del edge
template <typename TMPL>
    int del(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator& e);

//! @brief del face
template <typename TMPL>
    int del(mesh_tmpl<TMPL> &m, const typename TMPL::faces_t::const_iterator& f);

//! @brief get the number of edges of face f
template<typename TMPL>
size_t valence(const typename TMPL::faces_t::const_iterator &f);
//! @brief get the number of verts or edges of vert v
template<typename TMPL>
size_t valence(const typename TMPL::verts_t::const_iterator &v);

//! @brief set opposite and boundary edge
//! used after adding all faces, the result is a halfedge mesh
template <typename TMPL>
int set_opposite_and_boundary_edge(mesh_tmpl<TMPL> &m);

//! @brief copy a half-edge mesh to another
template <typename TMPL_A, typename TMPL_B>
int copy(const mesh_tmpl<TMPL_A>& ma, mesh_tmpl<TMPL_B>& mb);

//! @brief get the adjacent edges of face f
template <typename TMPL, typename CON>
void face_adj_edges(const typename TMPL::faces_t::const_iterator &f, CON& fe);
//! @brief get the adjacent verts of face f
template <typename TMPL, typename CON>
void face_adj_verts(const typename TMPL::faces_t::const_iterator &f, CON& fv);
//! @brief get the adjacent faces of face f
template <typename TMPL, typename CON>
void face_adj_faces(const typename TMPL::faces_t::const_iterator &f, CON& ff);
//! @brief get the adjacent faces of edge e
template <typename TMPL, typename CON>
void edge_adj_faces(const typename TMPL::edges_t::const_iterator &e, CON& ef);

//! @brief get the adjacent faces of vert v
template <typename TMPL, typename CON>
void vert_adj_faces(const typename TMPL::verts_t::const_iterator &v, CON& vf);

//! @brief get the adjacent verts of vert v
template <typename TMPL, typename CON>
void vert_adj_verts(const typename TMPL::verts_t::const_iterator &v, CON& vv);
//! @brief get the adjacent out edges of vert v
template <typename TMPL, typename CON>
void vert_adj_out_edges(const typename TMPL::verts_t::const_iterator &v, CON& ve);

//! @brief get the half edge of two vertices.
template <typename TMPL>
typename TMPL::edges_t::const_iterator
get_edge(const typename TMPL::verts_t::const_iterator & v1,
         const typename TMPL::verts_t::const_iterator & v2);

//! @brief gurantee that edge of boundary vertex is boundary edge
template <typename TMPL, typename CVI>
void
adjust_vert_edge(mesh_tmpl<TMPL> &m, const CVI & vert);

// @brief all edges in sec are border edge, and the region between
// sec[i*2] and sec[i*2+1] are connected by faces.
template <typename CVI, typename Con>
void sectors(const CVI &vi, Con &sec);

//! @brief add a new sector from in to out into sec., in and out are
//! outer edge of a sector, but may be not border edges.
template <typename TMPL, typename CEI, typename Con>
int add_face_into_sectors(mesh_tmpl<TMPL> &m, Con &sec, const CEI &in, const CEI &out);

//! @brief add a new face, hold topology traits after the operation
template <typename TMPL, typename RND_ITR>
typename TMPL::faces_t::const_iterator
add_face2(mesh_tmpl<TMPL> &m, const RND_ITR &vert_loop_beg, std::size_t n);

//! @brief is an given halfedge can be flipped
//! @param is_valid(i) = true
template <typename TMPL>
int try_edge_flip(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e);

//! @param can_edge_flip(m, i) == 0;
template <typename TMPL>
void edge_flip_by_rotate(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e);

//! @param can_edge_flip(m, i) == 0;
template <typename TMPL>
void edge_flip_by_del_add(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e);

//! @brief test whether an edge can be collapse
//! @param m for mesh, e for the edge to be collapsed
//! @return 0 for can, 1 for collapse may left some feature need to be
//! adjust by is_collapse_ok, 2 for can not collapse
template <typename TMPL>
int try_collapse(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e);


//! @brief implement a collapse operation on a given halfedge
//! @param m for mesh, e for the edge to be collapsed
//! @return the end vertex of the collapsed edge 
//! @warning {1. the operation will reserve the edges which have no faces
//! after collapse operation.
//! 2. The operation will try to keep traits of the mesh as its
//! original. i.e. if there is a hole in mesh before collapse, it will
//! be reserved}
template <typename TMPL>
const typename TMPL::verts_t::const_iterator&
collapse_edge(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e);

//! @brief post operation after a collapse to adjust some
//! non-manifold case e.g. multipal edges between two edges
//! @param the end vertex of the latest collapsed edge
//! @return 0 if operation
template <typename TMPL>
int is_collapse_ok(mesh_tmpl<TMPL> &m, const typename TMPL::verts_t::const_iterator &v);

//! @brief subdivide edge
template <typename TMPL>
const typename TMPL::verts_t::const_iterator
split_edge(mesh_tmpl<TMPL> &m, const typename TMPL::edges_t::const_iterator &e);

#include "operation.hh"

}}

#endif
