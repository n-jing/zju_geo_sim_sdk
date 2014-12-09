#include "halfedge_connectivity.h"
#include <iostream>

namespace gmesh {

std::vector<face_handle> halfedge_connectivity::get_adj_face(vert_handle vh) const {
  //assert(is_valid_handle(vh));
  halfedge_handle curr_heh = get_halfedge_handle(vh);
  assert(is_valid_handle(curr_heh));
  std::vector<face_handle> fhandles;
  if (!is_valid_handle(curr_heh)) return fhandles;

  face_handle fh = get_face_handle(curr_heh);
  while(!fh.is_valid()) { // boundary face_handle
    curr_heh = get_next_halfedge_handle(get_oppo_halfedge_handle(curr_heh));
    fh = get_face_handle(curr_heh);
  }
  do{
    if (fh.is_valid()) fhandles.push_back(fh);
    curr_heh = get_next_halfedge_handle(get_oppo_halfedge_handle(curr_heh));
    fh = get_face_handle(curr_heh);
  }while(curr_heh != get_halfedge_handle(vh) && is_valid_handle(curr_heh));
  return fhandles;
}

std::vector<face_handle> halfedge_connectivity::get_adj_face(edge_handle eh) const {
  halfedge_handle heh0 = get_halfedge_handle(eh, 0);
  halfedge_handle heh1 = get_halfedge_handle(eh, 1);
  assert(is_valid_handle(heh0) && is_valid_handle(heh1));
  face_handle fh0 = get_face_handle(heh0), fh1 = get_face_handle(heh1);
  std::vector<face_handle> fhandles;
  if(fh0.is_valid()) fhandles.push_back(fh0);
  if(fh1.is_valid()) fhandles.push_back(fh1);
  return fhandles;
}

std::vector<face_handle> halfedge_connectivity::get_adj_face(face_handle fh) const {
  std::vector<edge_handle> ehandles = get_edges(fh);
  std::vector<face_handle> fhandles;
  for(size_t i=0; i<ehandles.size(); ++i) {
    std::vector<face_handle> adj_fhs = get_adj_face(ehandles[i]);
    assert(adj_fhs.size() >0 && adj_fhs.size() <= 2);
    if(adj_fhs.size() == 2) {
      fhandles.push_back(adj_fhs[0] == fh ? adj_fhs[1] : adj_fhs[0]);
    }else if(adj_fhs.size() == 1 && adj_fhs[0] != fh)
      fhandles.push_back(adj_fhs[0]);
  }
  return fhandles;
}

std::vector<vert_handle> halfedge_connectivity::get_adj_vert(vert_handle vh) const {
  halfedge_handle curr_heh = get_halfedge_handle(vh);
  assert(is_valid_handle(curr_heh));
  std::vector<vert_handle> vhandles;
  do{
    vhandles.push_back(get_to_vert_handle(curr_heh));
    curr_heh = get_next_halfedge_handle(get_oppo_halfedge_handle(curr_heh));
  }while(curr_heh != get_halfedge_handle(vh) && !is_boundary(curr_heh));
  return vhandles;
}


std::vector<edge_handle> halfedge_connectivity::get_adj_edge(vert_handle vh) const {
  // if vh is a boundary vertex, the first edge should be a boundary edge
  halfedge_handle curr_heh = get_halfedge_handle(vh);
  assert(is_valid_handle(curr_heh));
  std::vector<edge_handle> ehandles;
  do{
    ehandles.push_back(get_edge_handle(curr_heh));
    curr_heh = get_next_halfedge_handle(get_oppo_halfedge_handle(curr_heh));
  }while(curr_heh != get_halfedge_handle(vh) && !is_boundary(curr_heh));
  return ehandles;
}

std::vector<vert_handle> halfedge_connectivity::get_verts(face_handle fh) const {
  halfedge_handle curr_heh = get_halfedge_handle(fh);
  assert(is_valid_handle(curr_heh));
  std::vector<vert_handle> vhandles;
  do{
    vhandles.push_back(get_from_vert_handle(curr_heh));
    curr_heh = get_next_halfedge_handle(curr_heh);
  }while(curr_heh != get_halfedge_handle(fh));
  return vhandles;
}

std::vector<edge_handle> halfedge_connectivity::get_edges(face_handle fh) const {
  halfedge_handle curr_heh = get_halfedge_handle(fh);
  assert(is_valid_handle(curr_heh));
  std::vector<edge_handle> ehandles;
  do{
    ehandles.push_back(get_edge_handle(curr_heh));
    curr_heh = get_next_halfedge_handle(curr_heh);
  }while(curr_heh != get_halfedge_handle(fh));
  return ehandles;
}

std::vector<vert_handle> halfedge_connectivity::get_verts(edge_handle eh) const {
  halfedge_handle heh0 = get_halfedge_handle(eh, 0);
  halfedge_handle heh1 = get_halfedge_handle(eh, 1);
  std::vector<vert_handle> vhandles;
  vhandles.push_back(get_vert_handle(heh0));
  vhandles.push_back(get_vert_handle(heh1));
  return vhandles;
}

halfedge_handle halfedge_connectivity::find_halfedge(vert_handle vh1, vert_handle vh2) const {
  assert(vh1.is_valid() && vh2.is_valid());
  halfedge_handle curr_heh = get_halfedge_handle(vh1);
  if(!is_valid_handle(curr_heh)) return halfedge_handle(-1);
  do{
    if(get_to_vert_handle(curr_heh) == vh2) return curr_heh;
    curr_heh = get_cw_rotated_halfedge_handle(curr_heh);
  }while(curr_heh != get_halfedge_handle(vh1) && curr_heh.is_valid());
  return halfedge_handle(-1); // cann't find a valid halfedge
}

edge_handle halfedge_connectivity::add_edge(vert_handle vh0, vert_handle vh1) {
  // create edge
  edge_handle eh = new_edge();
  halfedge_handle heh0 = new_halfedge();
  halfedge_handle heh1 = new_halfedge();
  set_halfedge_handle(eh, heh0, heh1);
  
  halfedge& he0 = get_halfedge(heh0);  
  he0.edge_handle_ = eh;
  he0.vert_handle_ = vh1;
  //he0.oppo_he_handle_ = heh1;

  halfedge& he1 = get_halfedge(heh1);
  he1.edge_handle_ = eh;
  he1.vert_handle_ = vh0;
  //he1.oppo_he_handle_ = heh0;
  return eh;
}

face_handle halfedge_connectivity::add_face(const std::vector<vert_handle>& vhandles) {
  size_t n = vhandles.size();
  // don't allow degenerated faces
  assert (n > 2);

  // Check sufficient working storage available
  if (edge_cache_.size() < n) {
    edge_cache_.resize(n);
    next_cache_.resize(6*n);
  }
  next_cache_count_ = 0;

  // test for topological errors
  for(size_t i=0, j=1; i<n; ++i, ++j, j%=n) {
    if( !is_boundary(vhandles[i])) {
      std::cerr<<  "PolyMeshT::add_face: complex vertex\n";
      std::cerr<< "The index of the vertex: " << vhandles[i].index()<< "\n";
      return face_handle(-1);
    }

    edge_cache_[i].he_handle = find_halfedge(vhandles[i], vhandles[j]);
    edge_cache_[i].is_new = !edge_cache_[i].he_handle.is_valid();
    edge_cache_[i].needs_adjust = false;

    if(!edge_cache_[i].is_new && !is_boundary(edge_cache_[i].he_handle)){
      std::cerr << "PolyMeshT::add_face: complex edge\n";
      return face_handle(-1);
    }
  }

  halfedge_handle inner_next, inner_prev;
  halfedge_handle outer_next, outer_prev;
  halfedge_handle bound_next, bound_prev;

  for(size_t i=0, j=1; i<vhandles.size(); ++i, ++j, j%=n) {    
    if(!edge_cache_[i].is_new && !edge_cache_[j].is_new) {
      // this two edges have been created
      inner_prev = edge_cache_[i].he_handle; // 
      inner_next = edge_cache_[j].he_handle;

      if(get_next_halfedge_handle(inner_prev) != inner_next) {
        outer_prev = get_oppo_halfedge_handle(inner_next);
        outer_next = get_oppo_halfedge_handle(inner_prev);
        bound_prev = outer_prev;
        do {
          bound_prev = get_oppo_halfedge_handle(get_next_halfedge_handle(bound_prev));
        } while(!is_boundary(bound_prev) || bound_prev == inner_prev);
        bound_next = get_next_halfedge_handle(bound_prev);
        assert(is_boundary(bound_prev) && is_boundary(bound_next));

        if(bound_next == inner_next) {
          std::cerr << "relink error\n";
          return face_handle(-1);
        }

        halfedge_handle p_start = get_next_halfedge_handle(inner_prev);
        // modified by fxz: next->prev
        halfedge_handle p_end = get_prev_halfedge_handle(inner_next);

        // store for cache
        next_cache_[next_cache_count_++] = std::make_pair(bound_prev, p_start);
        next_cache_[next_cache_count_++] = std::make_pair(p_end, bound_next);
        // next_cache_[next_cache_count_++] = std::make_pair(inner_prev, inner_next);	
      }
    }
  }

  // create new edges
  for(size_t i=0, j=1; i < vhandles.size(); ++i, ++j, j%=n) {
    if(edge_cache_[i].is_new) {
      edge_handle eh = add_edge(vhandles[i], vhandles[j]);
      edge_cache_[i].he_handle = get_halfedge_handle(eh, 0);
    }
  }

  // create new face
  face_handle fh(new_face());
  // set the last edge's first halfedge_handle as face's halfedge_handle
  set_halfedge_handle(fh, edge_cache_[n-1].he_handle);

  // setup halfedge
  for(size_t i=0, j=1; i<vhandles.size(); ++i, ++j, j%=n) {
    vert_handle vh = vhandles[j];
    inner_prev = edge_cache_[i].he_handle;
    inner_next = edge_cache_[j].he_handle;
    outer_prev = get_oppo_halfedge_handle(inner_next);
    outer_next = get_oppo_halfedge_handle(inner_prev);

    int cases = 0;
    if(edge_cache_[i].is_new) cases |= 0x01;
    if(edge_cache_[j].is_new) cases |= 0x02;

    switch(cases) {
      case 0: // both are old edges
        // a vertex's hh(halfedge_handle) should be a boundary hh
        // if(vh's hh is inner_next), but inner_next is not a boundary any more
        // so we need adjust vh's hh, and set a boundary hh to it if existing
        edge_cache_[j].needs_adjust = (get_halfedge_handle(vh) == inner_next);
        break;
      case 1: // prev is new, next is old
        bound_prev = get_prev_halfedge_handle(inner_next);
        next_cache_[next_cache_count_++] = std::make_pair(bound_prev, outer_next);
        set_halfedge_handle(vh, outer_next);
        break;
      case 2: // prev is old, next is new
        bound_next = get_next_halfedge_handle(inner_prev);
        next_cache_[next_cache_count_++] = std::make_pair(outer_prev, bound_next);
        set_halfedge_handle(vh, bound_next);
        break;	
      case 3: // both are new edges
        if( get_halfedge_handle(vh).is_valid()) {
          // this means there are other existing faces connecting vh, then
          // we need adjust the boundary halfedge connection
          bound_next = get_halfedge_handle(vh);
          bound_prev = get_prev_halfedge_handle(bound_next);
          next_cache_[next_cache_count_++] = std::make_pair(bound_prev, outer_next);
          next_cache_[next_cache_count_++] = std::make_pair(outer_prev, bound_next);	  
        }else { // this vertex has no other existing face
          set_halfedge_handle(vh, outer_next);
          next_cache_[next_cache_count_++] = std::make_pair(outer_prev, outer_next);
        }
        break;
      default:
        break;
    }
    next_cache_[next_cache_count_++] = std::make_pair(inner_prev, inner_next);
    set_face_handle(edge_cache_[i].he_handle, fh);
  }

  // process next/prev halfedge cache
  for(size_t i=0; i<next_cache_count_; ++i) {
    set_next_handle(next_cache_[i].first, next_cache_[i].second);
    set_prev_handle(next_cache_[i].second, next_cache_[i].first);
  }
  // adjust halfedge handle for boundary
  for(size_t i=0; i<n; ++i){
    if(edge_cache_[i].needs_adjust)
      adjust_vert_halfedge_handle(vhandles[i]);
  }

  return fh;
}

// add by fxz : only for trimesh
/*
void halfedge_connectivity::flip_edge(edge_handle eh)
{
  assert(!is_boundary(eh));
  
  halfedge_handle heha, hehb, hehan, hehbn;
  heha = get_halfedge_handle(eh, 0);
  hehb = get_halfedge_handle(eh, 1);

  face_handle fa, fb;
  fa = get_face_handle(heha);
  fb = get_face_handle(hehb);

  std::vector<vert_handle> vhl = get_verts(fa);
  assert(vhl.size() == 3);
  vhl.clear();
  vhl = get_verts(fb);
  assert(vhl.size() == 3);
  
  hehan = get_next_halfedge_handle(heha);
  hehbn = get_next_halfedge_handle(hehb);

  vert_handle va, vb, van, vbn;
  va = get_to_vert_handle(heha);
  vb = get_to_vert_handle(hehb);
  van = get_to_vert_handle(hehan);
  vbn = get_to_vert_handle(hehbn);

  delete_edge(eh);

  std::vector<vert_handle> vh;
  vh.resize(3);
  vh[0] = va; vh[1] = van; vh[2] = vbn;
  add_face(vh);
  vh[0] = vb; vh[1] = vbn; vh[2] = van;
  add_face(vh);
  
  return;
}
*/

bool halfedge_connectivity::is_flip_ok(edge_handle eh)
{
  if (is_boundary(eh)) return false;
  halfedge_handle h0, h1;
  h0 = get_halfedge_handle(eh, 0);
  h1 = get_halfedge_handle(eh, 1);
  
  vert_handle ah = get_to_vert_handle(get_next_halfedge_handle(h0));
  vert_handle bh = get_to_vert_handle(get_next_halfedge_handle(h1));

  if (ah == bh) return false;

  halfedge_handle heh = get_halfedge_handle(ah);
  do {
    if (get_to_vert_handle(heh) == bh)
      return false;
    heh = get_oppo_halfedge_handle(get_prev_halfedge_handle(heh));
  } while (heh != get_halfedge_handle(ah));
  
  return true;
}

void halfedge_connectivity::flip_edge(edge_handle eh)
{
  assert(is_flip_ok(eh));
  assert(!is_boundary(eh));
  
  halfedge_handle h0, h1, n0, n1, p0, p1;
  h0 = get_halfedge_handle(eh, 0);
  h1 = get_halfedge_handle(eh, 1);
  n0 = get_next_halfedge_handle(h0);
  n1 = get_next_halfedge_handle(h1);
  p0 = get_prev_halfedge_handle(h0);
  p1 = get_prev_halfedge_handle(h1);

  face_handle f0, f1;
  f0 = get_face_handle(h0);
  f1 = get_face_handle(h1);

  std::vector<vert_handle> vhl = get_verts(f0);
  assert(vhl.size() == 3);
  vhl.clear();
  vhl = get_verts(f1);
  assert(vhl.size() == 3);

  vert_handle v0, v1, v2, v3;
  
  v0 = get_to_vert_handle(h0);
  v1 = get_to_vert_handle(h1);
  v2 = get_to_vert_handle(n0);
  v3 = get_to_vert_handle(n1);

  // change some connectivity
  set_next_handle(n0, h0); set_prev_handle(h0, n0);
  set_next_handle(p1, n0); set_prev_handle(n0, p1);
  set_next_handle(h0, p1); set_prev_handle(p1, h0);
  set_next_handle(p0, n1); set_prev_handle(n1, p0);
  set_next_handle(n1, h1); set_prev_handle(h1, n1);
  set_next_handle(h1, p0); set_prev_handle(p0, h1);
  set_halfedge_handle(v0, n0);
  set_halfedge_handle(v1, n1);
  set_vert_handle(h0, v3);
  set_vert_handle(h1, v2);
  set_face_handle(p0, f1);
  set_face_handle(p1, f0);

  if (get_halfedge_handle(f0) == p0)
    set_halfedge_handle(f0, h0);
  if (get_halfedge_handle(f1) == p1)
    set_halfedge_handle(f1, h1);
  
  return;
}

bool halfedge_connectivity::is_collapse_ok(halfedge_handle heh)
{
  assert(is_valid_handle(heh));

  if (get_status(heh).is_deleted()) return false;
  
  vert_handle v0, v1;
  halfedge_handle oheh = get_oppo_halfedge_handle(heh);
  v0 = get_to_vert_handle(heh);
  v1 = get_to_vert_handle(oheh);
  
  halfedge_handle tmp = heh;
  std::vector<vert_handle> vv;
  vv.reserve(4);
  while (true) {
    tmp = get_oppo_halfedge_handle(get_prev_halfedge_handle(tmp));
    if (tmp == heh) break;
    vert_handle v = get_to_vert_handle(tmp), nv;
    
    // find the edge of v -> v0:
    halfedge_handle nh = get_halfedge_handle(v);
    do {
      nv = get_to_vert_handle(nh);
      nh = get_oppo_halfedge_handle(get_prev_halfedge_handle(nh));
    } while (nv!=v0 && nh!=get_halfedge_handle(v));
    
    if(nv==v0) vv.push_back(v);
  }
  
  if (vv.size() > 2) return false;

  if (vv.size() == 2) {
    halfedge_handle nh = get_halfedge_handle(vv[0]);
    do {
      if (get_to_vert_handle(nh) == vv[1]) return false;
      nh = get_next_halfedge_handle(get_oppo_halfedge_handle(nh));
    } while(nh != get_halfedge_handle(vv[0]));
  }
  
  return true;
}

void halfedge_connectivity::collapse_edge(halfedge_handle heh)
{
  assert(is_valid_handle(heh) && !get_status(heh).is_deleted());
  assert(is_collapse_ok(heh));
  
  edge_handle eh = get_edge_handle(heh);
  halfedge_handle n0, n1, p0, p1, h0, h1;
  vert_handle v0, v1;
  face_handle f0, f1;

  h0 = heh;
  h1 = get_oppo_halfedge_handle(heh);
  n0 = get_next_halfedge_handle(h0);
  n1 = get_next_halfedge_handle(h1);
  p0 = get_prev_halfedge_handle(h0);
  p1 = get_prev_halfedge_handle(h1);

  v0 = get_to_vert_handle(h0);
  v1 = get_to_vert_handle(h1);
  f0 = get_face_handle(h0);
  f1 = get_face_handle(h1);

  halfedge_handle tmp = p0;
  while(tmp != h1) {
    set_vert_handle(tmp, v0);
    tmp = get_prev_halfedge_handle(get_oppo_halfedge_handle(tmp));
  }
  set_next_handle(p0, n0); set_prev_handle(n0, p0);
  set_next_handle(p1, n1); set_prev_handle(n1, p1);
  if (connectivity::is_valid_handle(f0))
    set_halfedge_handle(f0, n0);
  if (connectivity::is_valid_handle(f1))
    set_halfedge_handle(f1, n1);
  set_halfedge_handle(v0, n0);

  if (connectivity::is_valid_handle(f0) &&  get_next_halfedge_handle(n0) == p0) {
    halfedge_handle on0, op0;
    on0 = get_oppo_halfedge_handle(n0);
    op0 = get_oppo_halfedge_handle(p0);
    face_handle fh_p = get_face_handle(op0);
    edge_handle eh_p = get_edge_handle(p0);
    halfedge_handle nop0, pop0;
    nop0 = get_next_halfedge_handle(op0);
    pop0 = get_prev_halfedge_handle(op0);
    set_next_handle(n0, nop0); set_prev_handle(nop0, n0);
    set_next_handle(pop0, n0); set_prev_handle(n0, pop0);
    set_halfedge_handle(v0, n0);
    set_face_handle(n0, fh_p);
    get_status(f0).set_deleted(true);

    vert_handle vh_n0 = get_to_vert_handle(n0);
    set_halfedge_handle(vh_n0, on0);
      
    if (connectivity::is_valid_handle(fh_p)) {
      set_halfedge_handle(fh_p, n0);
    }

    get_status(p0).set_deleted(true);
    get_status(op0).set_deleted(true);
    get_status(eh_p).set_deleted(true);
      
    if (nop0 == on0) {
      get_status(vh_n0).set_deleted(true);
      get_status(n0).set_deleted(true);
      get_status(on0).set_deleted(true);
      edge_handle eh_n0 = get_edge_handle(n0);
      get_status(eh_n0).set_deleted(true);
      halfedge_handle non0 = get_next_halfedge_handle(on0);
      if (non0 == n0) {
        set_halfedge_handle(v0, halfedge_handle(-1));
      } else {
        set_next_handle(pop0, non0); set_prev_handle(non0, pop0);
        set_halfedge_handle(v0, non0);
      }
    }
  }

  if (connectivity::is_valid_handle(f1) && get_next_halfedge_handle(n1) == p1) {
    halfedge_handle op1, on1;
    op1 = get_oppo_halfedge_handle(p1);
    on1 = get_oppo_halfedge_handle(n1);
    face_handle fh_p1 = get_face_handle(op1);
    edge_handle eh_p1 = get_edge_handle(p1);
    halfedge_handle nop1, pop1;
    nop1 = get_next_halfedge_handle(op1);
    pop1 = get_prev_halfedge_handle(op1);
    set_next_handle(n1, nop1); set_prev_handle(nop1, n1);
    set_next_handle(pop1, n1); set_prev_handle(n1, pop1);
    set_halfedge_handle(v0, n1);
    set_face_handle(n1, fh_p1);
    get_status(f1).set_deleted(true);

    vert_handle vh_n1 = get_to_vert_handle(n1);
    set_halfedge_handle(vh_n1, on1);
    if (connectivity::is_valid_handle(fh_p1)) {
      set_halfedge_handle(fh_p1, n1);
    }

    get_status(p1).set_deleted(true);
    get_status(op1).set_deleted(true);
    get_status(eh_p1).set_deleted(true);

    if (nop1 == on1) {
      get_status(vh_n1).set_deleted(true);
      get_status(n1).set_deleted(true);
      get_status(on1).set_deleted(true);
      edge_handle eh_n1 = get_edge_handle(n1);
      get_status(eh_n1).set_deleted(true);
      halfedge_handle non1 = get_next_halfedge_handle(on1);
      if (non1 == n1) {
        set_halfedge_handle(v0, halfedge_handle(-1));
      } else {
        set_next_handle(pop1, non1); set_prev_handle(non1, pop1);
        set_halfedge_handle(v0, non1);
      }
    }
  }

  get_status(eh).set_deleted(true);
  get_status(h0).set_deleted(true);
  get_status(h1).set_deleted(true);
  get_status(v1).set_deleted(true);

  if (!is_valid_handle(get_halfedge_handle(v0)))
    get_status(v0).set_deleted(true);
  
  return;
}

void halfedge_connectivity::adjust_vert_halfedge_handle(vert_handle vh) {
  halfedge_handle he = get_halfedge_handle(vh);
  do{
    if(is_boundary(he)) {
      set_halfedge_handle(vh, he);
      break;
    }
    he = get_cw_rotated_halfedge_handle(he);
  }while( he != get_halfedge_handle(vh));
}

void halfedge_connectivity::delete_vert(vert_handle vh) {
  std::vector<face_handle> fhs; /// store adj faces
  /*fhs.reserve(8);
  for(vf_iter vf_it(vf_iter(*this, vh)); vf_it; ++vf_it) {
    fhs.push_back(vf_it.handle());
  }*/
  fhs = get_adj_face(vh);

  for(std::vector<face_handle>::iterator fh_it = fhs.begin();
      fh_it != fhs.end(); ++fh_it) {
    delete_face(*fh_it, true);
  }

  get_status(vh).set_deleted(true);
}

void halfedge_connectivity::delete_edge(edge_handle eh, bool is_del_iso_vert) {
  std::vector<face_handle> fhandles = get_adj_face(eh);
  for(size_t i=0; i<fhandles.size(); ++i)
    delete_face(fhandles[i], is_del_iso_vert);
}

void halfedge_connectivity::delete_face(face_handle fh, bool is_del_iso_vert) {
  //assert(fh.is_valid() && !(get_status(fh).is_deleted()));
  // modify by fxz : used for test
  assert(fh.is_valid());
  assert(!(get_status(fh).is_deleted()));
  
  get_status(fh).set_deleted(true);
  
  std::vector<edge_handle> ehandles; ehandles.reserve(3);
  std::vector<vert_handle> vhandles; vhandles.reserve(3);

  /// for all halfedges of face fh, do:
  ///   1) invalidate face handle
  ///   2) collect all boundary halfedges, set them deleted
  ///   3) store vertex handles
  for(fh_iter fh_it(fh_iter(*this, fh)); fh_it; ++fh_it) {
    halfedge_handle hh = fh_it.handle();    
    set_boundary(hh);
    if(is_boundary(get_oppo_halfedge_handle(hh)))
      ehandles.push_back(get_edge_handle(hh));
    vhandles.push_back(get_to_vert_handle(hh));
  }  

  /// delete all collected edges
  for(size_t i=0; i<ehandles.size(); ++i) {
    edge_handle eh = ehandles[i];
    halfedge_handle h0 = get_halfedge_handle(eh, 0);
    halfedge_handle n0 = get_next_halfedge_handle(h0);
    halfedge_handle p0 = get_prev_halfedge_handle(h0);
    vert_handle v0 = get_to_vert_handle(h0);
    
    halfedge_handle h1 = get_halfedge_handle(eh, 1);
    halfedge_handle n1 = get_next_halfedge_handle(h1);
    halfedge_handle p1 = get_prev_halfedge_handle(h1);
    vert_handle v1 = get_to_vert_handle(h1);

    // add by fxz
    set_prev_handle(n1, p0);
    set_prev_handle(n0, p1);

    set_next_handle(p0, n1);
    set_next_handle(p1, n0);

    /// make edge/halfedge deleted
    get_status(eh).set_deleted(true);
    get_status(h0).set_deleted(true);
    get_status(h1).set_deleted(true);

    /// update v0, v1
    if(get_halfedge_handle(v0) == h1) {
      if(n0 == h1) {
        if (is_del_iso_vert)  // judge whether to delete the isolated vert
          get_status(v0).set_deleted(true);
        set_halfedge_handle(v0, halfedge_handle(-1));
      }else
        set_halfedge_handle(v0, n0);
    }
    if(get_halfedge_handle(v1) == h0) {
      if(n1 == h0) {
        if (is_del_iso_vert)  // judge whether to delete the isolated vert
          get_status(v1).set_deleted(true);
        set_halfedge_handle(v1, halfedge_handle(-1));
      }else
        set_halfedge_handle(v1, n1);
    }
  }

  for(size_t i=0; i<vhandles.size(); ++i)
    if (get_halfedge_handle(vhandles[i]).is_valid())
      adjust_vert_halfedge_handle(vhandles[i]);
}

void halfedge_connectivity::garbage_collection() {
  size_t nv = pm_.n_verts(), ne = pm_.n_edges();
  size_t nf = pm_.n_faces(), nh = pm_.n_halfedges();
  std::vector<vert_handle> vh_map; vh_map.reserve(nv);
  std::vector<edge_handle> eh_map; eh_map.reserve(ne);
  std::vector<face_handle> fh_map; fh_map.reserve(nf);
  std::vector<halfedge_handle> hh_map; hh_map.reserve(nh);
  
  for(int i=0; i<nv; ++i) vh_map.push_back(vert_handle(i));
  for(int i=0; i<ne; ++i) eh_map.push_back(edge_handle(i));
  for(int i=0; i<nf; ++i) fh_map.push_back(face_handle(i));
  for(int i=0; i<nh; ++i) hh_map.push_back(halfedge_handle(i));

  /// remove deleted verts
  if(nv > 0) {
    int i0 = 0, i1 = nv -1;
    while(1) {
      while(!get_status(vert_handle(i0)).is_deleted() && i0 < i1) ++i0;
      // modify by fxz
      //while(!get_status(vert_handle(i1)).is_deleted() && i0 < i1) --i1;
      while(get_status(vert_handle(i1)).is_deleted() && i0 < i1) --i1;
      //if(i0 > i1) break;
      if (i0 >= i1) break;

      std::swap(verts_[i0], verts_[i1]);
      std::swap(vh_map[i0], vh_map[i1]);
      pm_.vprops_swap(i0, i1);
    }

    verts_.resize(get_status(vert_handle(i0)).is_deleted() ? i0 : i0 + 1);
    pm_.vprops_resize(verts_.size());    
  }
  /// remove deleted edges
  if(ne > 0) {
    int i0 = 0, i1 = ne -1;
    while(1) {
      while(!get_status(edge_handle(i0)).is_deleted() && i0 < i1) ++i0;
      //while(!get_status(edge_handle(i1)).is_deleted() && i0 < i1) --i1;
      //if(i0 > i1) break;
      // modify by fxz
      while(get_status(edge_handle(i1)).is_deleted() && i0 < i1) --i1;
      if (i0 >= i1) break;

      std::swap(edges_[i0], edges_[i1]);
      std::swap(eh_map[i0], eh_map[i1]);
      pm_.eprops_swap(i0, i1);      
    }

    edges_.resize(get_status(edge_handle(i0)).is_deleted() ? i0 : i0 + 1);

    // add by fxz
    //pm_.eprops_clear();

    pm_.eprops_resize(edges_.size());
  }
  /// remove deleted faces
  if(nf > 0) {
    int i0 = 0, i1 = nf -1;
    while(1) {
      while(!get_status(face_handle(i0)).is_deleted() && i0 < i1) ++i0;
      //while(!get_status(face_handle(i1)).is_deleted() && i0 < i1) --i1;
      //if(i0 > i1) break;
      // modify by fxz
      while (get_status(face_handle(i1)).is_deleted() && i0 < i1) --i1;
      if (i0 >= i1) break;

      std::swap(faces_[i0], faces_[i1]);
      std::swap(fh_map[i0], fh_map[i1]);
      pm_.fprops_swap(i0, i1);
    }

    faces_.resize(get_status(face_handle(i0)).is_deleted() ? i0 : i0 + 1);

    // add by fxz
    //pm_.fprops_clear();

    pm_.fprops_resize(faces_.size());
  }
  /// remove delete halfedges
  if(nh > 0) {
    int i0 = 0, i1 = nh -1;
    while(1) {
      while(!get_status(halfedge_handle(i0)).is_deleted() && i0 < i1) ++i0;
      //while(!get_status(halfedge_handle(i1)).is_deleted() && i0 < i1) --i1;
      //if(i0 > i1) break;
      // modify by fxz
      while (get_status(halfedge_handle(i1)).is_deleted() && i0 < i1) --i1;
      if (i0 >= i1) break;

      std::swap(halfedges_[i0], halfedges_[i1]);
      std::swap(hh_map[i0], hh_map[i1]);
      pm_.hprops_swap(i0, i1);
    }

    halfedges_.resize(get_status(halfedge_handle(i0)).is_deleted() ? i0 : i0 + 1);

    // add by fxz
    //pm_.hprops_clear();

    pm_.hprops_resize(halfedges_.size());
  }

  /// update vertex handle
  for(size_t i=0; i<verts_.size(); ++i) {
    vert_handle vh = connectivity::get_handle(verts_[i]);
    halfedge_handle heh = get_halfedge_handle(vh);
    if (heh.is_valid())
      set_halfedge_handle(vh, hh_map[heh.index()]);
  }
  /// update edge/halfedge handle
  for(size_t i=0; i<edges_.size(); ++i) {
    /// update edge handle first
    edge_handle eh = connectivity::get_handle(edges_[i]);
    set_halfedge_handle(eh, hh_map[get_halfedge_handle(eh, 0).index()],
                        hh_map[get_halfedge_handle(eh, 1).index()]);
    /// update halfedge handle
    halfedge_handle hh = get_halfedge_handle(eh, 0);
    set_vert_handle(hh, vh_map[get_to_vert_handle(hh).index()]);
    // add by fxz
    set_edge_handle(hh, eh);
    hh = get_halfedge_handle(eh, 1);
    set_vert_handle(hh, vh_map[get_to_vert_handle(hh).index()]);
    // add by fxz
    set_edge_handle(hh, eh);
  }
  /// update halfedge connectivity
  for(size_t i=0; i<edges_.size(); ++i) {
    edge_handle eh = connectivity::get_handle(edges_[i]);
    halfedge_handle hh = get_halfedge_handle(eh, 0);
    set_next_handle(hh, hh_map[get_next_halfedge_handle(hh).index()]);
    set_prev_handle(hh, hh_map[get_prev_halfedge_handle(hh).index()]);
    


    // add by fxz    
    //set_oppo_handle(hh, hh_map[get_oppo_halfedge_handle(hh).index()]);

    if(!is_boundary(hh)) {
      set_face_handle(hh, fh_map[get_face_handle(hh).index()]);
    }
    hh = get_halfedge_handle(eh, 1);
    set_next_handle(hh, hh_map[get_next_halfedge_handle(hh).index()]);
    set_prev_handle(hh, hh_map[get_prev_halfedge_handle(hh).index()]);

    // add by fxz    
    //set_oppo_handle(hh, hh_map[get_oppo_halfedge_handle(hh).index()]);

    if(!is_boundary(hh)) {
      set_face_handle(hh, fh_map[get_face_handle(hh).index()]);
    }
  }
  /// update face handle
  for(size_t i=0; i<faces_.size(); ++i) {
    face_handle fh = connectivity::get_handle(faces_[i]);
    set_halfedge_handle(fh, hh_map[get_halfedge_handle(fh).index()]);
  }
}
  
}
