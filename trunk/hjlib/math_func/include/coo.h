#ifndef HJ_MATH_FUNC_COO_H_
#define HJ_MATH_FUNC_COO_H_

#include <vector>
#include <algorithm>
#include <iostream>

#include "config.h"

namespace hj { namespace math_func {

//! @brief

inline size_t powi(size_t base, size_t pow) {
  size_t rtn = 1;
  for(size_t i = 0; i < pow; ++i)
    rtn *= base;
  return rtn;
}

//! c is dim major. TODO: put the requirement to function.
template <typename T>
void pack_coo(size_t num, size_t dim, const T *c, std::vector<std::vector<T> > &pc)
{
  using namespace std;
  pc.resize(num); {
    size_t ni;
#if HJ_FUNC_USE_OMP
#pragma omp parallel for private(ni)
#endif
    for(ni = 0; ni < num; ++ni)
      pc[ni].insert(pc[ni].begin(), c+ni*dim, c+ni*dim+dim);
    sort(pc.begin(), pc.end());
    typename vector<vector<T> >::iterator pos = unique(pc.begin(), pc.end());
    pc.erase(pos, pc.end());
  }
}

template <typename T> class coo_pat_coo;
template <typename T> class coo_pat_dense;
template <typename T> class coo_pat_csc;
class coo_l2g;
template <typename INT> class coo_pat;

// dim major
template <typename INT>
class coo_set
{
public:
  coo_set()
    :cur_pos_(0), format_(0) {
  }
  coo_set(size_t dim, size_t nf, size_t nx, size_t nnz)
    :cur_pos_(0) {
    if(nnz != -1) {
      coos_.resize(nnz*dim);
      format_ = 'S';
    }
    else {
      nnz = nf * powi(nx, dim-1);
      format_ = 'D';
    }
    dim_ = dim;
    nf_ = nf;
    nx_ = nx;
  }
  size_t nf(void) const { return nf_; }
  size_t nx(void) const { return nx_; }
  size_t dim(void) const { return dim_; }
  size_t coo_num(void) const { return coos_.size()/dim_; }
  bool is_dense(void) const { return format_ == 'D'; }

  coo_pat<INT> *get_coo_pat(void) const {
    switch(format_) {
    case 'D': return new coo_pat_dense<INT>(dim(), nf(), nx());
    case 'S': {
      std::vector<std::vector<INT> > pc;
      pack_coo(coo_num(), dim(), &coos_[0], pc);
      coo_pat_dense<INT> cpd(dim()-1, nf(), nx());
      if(dim() > 1 && cpd.nnz() < pc.size())
        return new coo_pat_csc<INT>(dim(), nf(), nx(), pc);
      //      std::cout << "use coo: " << pc.size() << " " << cpd.nnz() << std::endl;
      return new coo_pat_coo<INT>(dim(), nf(), nx(), pc);
    }
    }
    return 0;
  }
private:
  friend class coo_l2g;
  inline void add(const INT*c) {
    assert(!is_dense());
    std::copy(c, c+dim(), &coos_[cur_pos_]);
    cur_pos_ += dim();
  }
  friend class coo_pat_coo<INT>;
  friend class coo_pat_csc<INT>;
  std::vector<INT> coos_;
  size_t cur_pos_, dim_, nf_, nx_;
  char format_;
};

template <typename INT> class coo_pat_dense;
template <typename INT> class coo_pat_coo;

template <typename INT>
class coo_pat
{
public:
  coo_pat(size_t dim, size_t nf, size_t nx, char format)
    :dim_(dim), nf_(nf), nx_(nx), format_(format), d_(0), s_(0), p_(0) {
    switch(format_) {
    case 'D': d_ = static_cast<const coo_pat_dense<INT> *>(this); break;
    case 'S': s_ = static_cast<const coo_pat_coo<INT> *>(this); break;
    case 'P': p_ = static_cast<const coo_pat_csc<INT> *>(this); break;
    }
  }

  inline size_t is_dense(void) const { return format_ == 'D'; }
  inline size_t nnz(void) const { return nnz_; }
  inline size_t dim(void) const { return dim_; }
  inline size_t nf(void) const { return nf_; }
  inline size_t nx(void) const { return nx_; }
  inline size_t operator()(const INT *coo) const {
    switch(format_) {
    case 'D': return (*d_)(coo);
    case 'S': return (*s_)(coo);
    case 'P': return (*p_)(coo);
    }
    return -1;
  }
  inline INT *operator()(size_t nzi, INT *coo) const {
    switch(format_) {
    case 'D': return (*d_)(nzi, coo);
    case 'S': return (*s_)(nzi, coo);
    case 'P': return (*p_)(nzi, coo);
    }
    return 0;
  }

  virtual ~coo_pat(){}
protected:
  const size_t dim_, nf_, nx_;
  size_t nnz_;
  const char format_;
  const coo_pat_dense<INT> *d_;
  const coo_pat_coo<INT> *s_;
  const coo_pat_csc<INT> *p_;
};

template <typename INT>
class coo_pat_dense : public coo_pat<INT>
{
public:
  using coo_pat<INT>::dim;
  using coo_pat<INT>::nx;
  using coo_pat<INT>::nnz_;
  coo_pat_dense(size_t dim, size_t nf, size_t nx)
    :coo_pat<INT>(dim, nf, nx, 'D'){
    nnz_ = nf * powi(nx, dim-1);
  }
  inline size_t operator()(const INT *coo) const {
    size_t ad = coo[0];
    for(size_t d = 1; d < dim(); ++d) {
      ad *= nx();
      ad += coo[d];
    }
    return ad;
  }
  inline INT *operator()(size_t nzi, INT *coo) const {
    for(size_t d = dim()-1; d > 0; --d) {
      coo[d] = nzi%nx();
      nzi = (nzi-coo[d])/nx();
    }
    coo[0] = nzi;
    return coo;
  }
};

template <typename INT>
class coo_pat_coo : public coo_pat<INT>
{
public:
  using coo_pat<INT>::dim;
  using coo_pat<INT>::nnz_;
  using coo_pat<INT>::nnz;
  coo_pat_coo(size_t dim, size_t nf, size_t nx, const std::vector<std::vector<INT> > &pc)
    :coo_pat<INT>(dim, nf, nx, 'S') {
    coos_.resize(dim*pc.size());
    for(size_t i = 0; i < pc.size(); ++i) {
      for(size_t d = 0; d < dim; ++d)
        coos_[i+d*pc.size()] = pc[i][d];
    }
    nnz_ = pc.size();
  }
  inline size_t operator()(const INT *coo) const {
    size_t beg = 0, end = nnz();
    const INT *ptr = &coos_[0];
    for(size_t i = 0; i < dim()-1; ++i, ptr+=nnz()) {
      const std::pair<const INT *, const INT *>
        r = std::equal_range(ptr+beg, ptr+end, coo[i]);
      beg = r.first - ptr;
      end = r.second - ptr;
    }
    beg = std::lower_bound(ptr+beg, ptr+end, coo[dim()-1])-ptr;
    if(beg == end || ptr[beg] != coo[dim()-1]) return -1;
    return beg;
  }
  inline INT *operator()(size_t nzi, INT *coo) const {
    assert(nzi < nnz());
    for(size_t i = 0; i < dim(); ++i, nzi += nnz())
      coo[i] = coos_[nzi];
    return coo;
  }
private:
  std::vector<INT> coos_;
};

template <typename INT>
class coo_pat_csc : public coo_pat<INT>
{
public:
  using coo_pat<INT>::dim;
  using coo_pat<INT>::nnz_;
  using coo_pat<INT>::nnz;
  coo_pat_csc(size_t dim, size_t nf, size_t nx, const std::vector<std::vector<INT> > &pc)
    :coo_pat<INT>(dim, nf, nx, 'P'),
     ptr_addr_(dim-1, nf, nx) {
    assert(dim > 1);
    nnz_ = pc.size();
    // the last is coo, and the prev are ptrs
    ptr_.resize(ptr_addr_.nnz()+1, 0);
    coos_.resize(pc.size());

    for(size_t pi = 1, nzi = 0; pi < ptr_.size(); ++pi) {
      ptr_[pi] = ptr_[pi-1];
      for(; nzi < pc.size(); ++nzi) {

        coos_[nzi] = pc[nzi][dim-1];

        const size_t pa = ptr_addr_(&pc[nzi][0]);
        assert(pa >= pi-1);
        if(pa == pi-1)
          ++ptr_[pi];
        else
          break;
      }
    }

#if 0
    std::vector<INT> c(dim());
    size_t nzi, nzi2;
    for(nzi = 0; nzi < nnz(); ++nzi) {
      (*this)(nzi, &c[0]);
      if(c != pc[nzi])
        std::cout << "error." << nzi << std::endl;
      nzi2 = (*this)(&c[0]);
      if(nzi2 != nzi)
        std::cout << "error." << nzi2 << std::endl;        
    }

#endif
  }
  inline size_t operator()(const INT *coo) const {
    const size_t pa = ptr_addr_(coo);
    size_t beg = ptr_[pa], end = ptr_[pa+1];
    const INT *idx = &coos_[0];
    size_t off = std::lower_bound(idx+beg, idx+end, coo[dim()-1])-idx;
    if(off == end || idx[off] != coo[dim()-1]) return -1;
    return off;
  }
  inline INT *operator()(size_t nzi, INT *coo) const {
    assert(nzi < nnz());
    coo[dim()-1] = coos_[nzi];
    const size_t pa = std::upper_bound(ptr_.begin(), ptr_.end(), nzi)-ptr_.begin()-1;
    assert(pa != -1);
    ptr_addr_(pa, coo);
    return coo;
  }
private:
  coo_pat_dense<INT> ptr_addr_;
  std::vector<INT> ptr_;
  std::vector<INT> coos_;
};

//! return value for a coodindates with coordinate transform
class coo_map {
public:
  virtual ~coo_map(){};
  coo_map(size_t f_base = 0)
    :f_base_(f_base) {
  }
  inline const size_t f_base(void) const { return f_base_; }
  inline size_t &f_base(void) { return f_base_; }
private:
  size_t f_base_;
};

template <typename VAL_TYPE, typename INT_TYPE>
class coo2val_t : public coo_map {
public:
  typedef VAL_TYPE val_type;
  typedef INT_TYPE int_type;
  coo2val_t(const coo_pat<int_type> &cp, val_type *val, size_t f_base = 0)
    :coo_map(f_base), cp_(cp), val_(val) {
  }
  //! NOTICE: coo will be changed
  inline val_type &operator[](int_type *coo) const {
    return val_[addr(coo)];
  }
private:
  inline size_t addr(int_type *coo) const {
    coo[0] += f_base();
    const size_t ad = cp_(coo);
    assert(ad != -1);
    return ad;
  }
  const coo_pat<int_type> &cp_;
  val_type *val_;
};

template <typename VAL_TYPE, typename INT_TYPE>
coo2val_t<VAL_TYPE, INT_TYPE>
coo2val(const coo_pat<INT_TYPE> &cp, VAL_TYPE *val, size_t f_base = 0) {
  return coo2val_t<VAL_TYPE, INT_TYPE>(cp, val, f_base);
}

class coo_l2g : public coo_map {
public:
  coo_l2g(size_t f_base = 0)
    :coo_map(f_base) {
  }
  //! NOTICE: c will be changed
  template <typename INT_TYPE>
  inline void add(coo_set<INT_TYPE> &cs, INT_TYPE *c) const {
    c[0] += f_base();
    cs.add(c);
  }
};

template <typename INT>
void coo2csc(const coo_pat<INT> &coo,  INT *ptr, INT *idx,
             size_t beg, size_t nnz, size_t dim_base)
{
  assert(coo.dim() >= 2+dim_base);
  ptr[0] = 0;
  std::vector<INT> c(coo.dim());
  for(size_t ti = 1, ni = beg; ti != -1;) {
    ptr[ti] = ptr[ti-1];
    for(; ni < nnz; ++ni) {
      if(coo(ni, &c[0])[dim_base] == ti-1)
        ++ptr[ti];
      else {
        ++ti;
        break;
      }
      idx[ni] = c[dim_base+1];
    }
    if(ni == nnz)
      ti = -1;
  }
}

// the last two dimension
template <typename INT>
void coo2csc(const coo_pat<INT> &coo,  INT *ptr, INT *idx, INT *leading = 0)
{
  assert(!coo.is_dense());
  if(leading == 0) {
    coo2csc(coo, ptr, idx, 0, coo.nnz(), 0);
    assert(ptr[coo.nf()] == coo.nnz());
  }
  else {
    assert(leading);
    if(coo.nf() == 1 && coo.dim() == 3)
      coo2csc(coo, ptr, idx, 0, coo.nnz(), 1);
    else
      assert(0); // may need turn coo into an iteratable object
  }
}

}}

#endif
