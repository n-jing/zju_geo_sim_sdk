#ifndef HJ_HALF_EDGE_CONTAINER_H_
#define HJ_HALF_EDGE_CONTAINER_H_

//! @file container.h
//! @brief basic container for half_edge
//! @author Jin Huang
//! @date 201406-

#include <list>
#include <vector>
#include <memory>
#include <cassert>

namespace hj { namespace half_edge {

//! @brief kinds of set, but can be unordered.  The element in it
//! should have constant iterator, thus std::vector::iterator does not
//! work.
template <typename T>
class container_concept
{
public:
  typedef T value_type; //! value type of container
  class iterator{};
  class const_iterator{};

  std::size_t size(void) const;
  const_iterator add(const value_type &v) { return const_iterator(); }
  void del(const const_iterator &i) { }
};

//! @brief wrap template argument for container
template <typename TMPL>
class cont_tmpl
{
public:
  typedef TMPL real_t;
	real_t &operator()(void) {return *static_cast<real_t *>(this);}
	const real_t &operator()(void) const {return *static_cast<const real_t *>(this);}
};

//! @brief wrap template argument for iterator
template <typename TMPL>
class iterator_tmpl
{
public:
 typedef TMPL real_t;
	real_t &operator()(void) {return *static_cast<real_t *>(this);}
	const real_t &operator()(void) const {return *static_cast<const real_t *>(this);}
  operator bool() const { return !!(*this); }
};

//! @brief return itr == null, makes itr look like a pointer
template <typename TMPL>
bool operator!(const iterator_tmpl<TMPL> &itr) { return itr() == TMPL(); }

/// begin of the real implementation

//! @brief container based on std::list
template <typename T>
class std_list : public std::list<T>, public cont_tmpl<std_list<T> >
{
public:
  typedef T value_type;

  //! @brief const_iterator
  class const_iterator : public std::list<T>::const_iterator,
                         public iterator_tmpl<const_iterator>
  {
  public:
    const_iterator(){}
    const_iterator(const typename std::list<T>::const_iterator &i):std::list<T>::const_iterator(i) {}
    const_iterator(const typename std::list<T>::iterator &i):std::list<T>::const_iterator(i) {}
  };

  //! @brief iterator
  class iterator : public std::list<T>::iterator,
                   public iterator_tmpl<iterator>
  {
  public:
    iterator(){}
    iterator(const typename std::list<T>::iterator &i):std::list<T>::iterator(i) {}
  };

  const_iterator add(const value_type &v) {
    base_t::push_back(v);
    const_iterator last(base_t::end());
    --last;
    return last;
  }

  void del(const const_iterator &i) {
    base_t::erase(cast_iterator<iterator>(*this, i));
  }

private:
  typedef typename std::list<T> base_t;

  //! @brief: bruteforce way to hack current compiler's complain.
  template <typename ITR, typename C, typename C_ITR>
  ITR cast_iterator(C &c, C_ITR i) {
    /// Safer but slow way
    // ITR it = c.begin();
    // std::advance (it, std::distance<C_ITR>(it, i));  
    // return it;
    return *reinterpret_cast<ITR *>(const_cast<void *>(reinterpret_cast<const void *>(&i)));
  }
};

//! @brief use pointer to indicate deleted or not.  For large internal
//! property, it is more memory efficient.
template <typename T>
class vector_by_ptr
{
public:
  T &operator[](size_t i) { assert(!is_deleted(i)); return *std_vec_[i]; }
  const T &operator[](size_t i) const { assert(!is_deleted(i)); return *std_vec_[i]; }
  void push_back(const T &t) {
    std_vec_.push_back(std::shared_ptr<T>(new T(t)));
  }
  void del(size_t i) { std_vec_[i].reset(); }
  bool is_deleted(size_t i) const { return !std_vec_[i].get(); }
  size_t size() const { return std_vec_.size(); }
private:
  std::vector<std::shared_ptr<T> > std_vec_;
};

//! @brief use status to indicate deleted or not.  For small internal
//! property, it is more time efficient.
//! @note which is the best? bitset or vector<int>, vector<char>?
template <typename T>
class vector_by_bit
{
public:
  T &operator[](size_t i) { assert(!is_deleted(i)); return std_vec_[i]; }
  const T &operator[](size_t i) const { assert(!is_deleted(i)); return std_vec_[i]; }
  void push_back(const T &t) {
    std_vec_.push_back(t);
    is_deleted_.push_back(false);
  }
  void del(size_t i) { is_deleted_[i] = true; }
  bool is_deleted(size_t i) const { return is_deleted_[i]; }
  size_t size() const { return std_vec_.size(); }
private:
  std::vector<T> std_vec_;
  std::vector<bool> is_deleted_;
};

//! @brief container based on std::vector. Be careful the reallocation
//! of std::vector makes the iterator invalid.
//! @todo can have different methods: shared_ptr<T>, pair<bool, T>,
//! etc.
template <typename T, template <typename T> class strategy>
class std_vector_base : public cont_tmpl<std_vector_base<T, strategy> >
{
private:
  typedef strategy<T> vector_t;
public:
  typedef T value_type;

  std_vector_base():size_(0){}

  size_t size() const { return size_; }

  //! @brief base for const and non-const iterator
  template <typename OWNER>
  class iterator_base : public std::iterator<std::bidirectional_iterator_tag, value_type>
  {
  public:
    iterator_base():addr_(-1), this_(0) {}
    iterator_base(std::size_t addr, OWNER &owner):addr_(addr), this_(&owner) {}

    const value_type *operator->() const { return &this_->std_vec_[addr_]; }
    const value_type &operator*() const { return this_->std_vec_[addr_]; }

    iterator_base<OWNER> operator++() {
      for(++addr_; addr_ < this_->std_vec_.size() && this_->std_vec_.is_deleted(addr_); ++addr_);
      return *this;
    }
    iterator_base<OWNER> operator--() {
      for(--addr_; this_->std_vec_.is_deleted(addr_) && addr_ > 0; --addr_);
      return *this;
    }
    template <typename O>
    bool operator == (const iterator_base<O> &i) const { return addr_ == i.addr_ && this_ == i.this_; }
    template <typename O>
    bool operator != (const iterator_base<O> &i) const { return !(*this == i); }
  protected:
    friend class std_vector_base<T, strategy>;
    std::size_t addr_;
    OWNER *this_;
  };

  //! @brief const_iterator
  class const_iterator : public iterator_base<const std_vector_base<T, strategy> >,
                         public iterator_tmpl<const_iterator>
  {
  public:
    const_iterator() {}
    const_iterator(std::size_t addr, const std_vector_base<T, strategy> &owner)
      :iterator_base<const std_vector_base<T, strategy> >(addr, owner) {}
  private:
    using iterator_base<const std_vector_base<T, strategy> >::this_;
    using iterator_base<const std_vector_base<T, strategy> >::addr_;
  };

  //! @brief iterator
  class iterator : public iterator_base<std_vector_base<T, strategy> >,
                   public iterator_tmpl<iterator>
  {
  public:
    iterator(){}
    iterator(std::size_t addr, std_vector_base<T, strategy> &owner)
      :iterator_base<std_vector_base<T, strategy> >(addr, owner) {}
    value_type *operator->() { return &this_->std_vec_[addr_]; }
    value_type &operator*() { return this_->std_vec_[addr_]; }
    operator const_iterator() const { return const_iterator(addr_, *this_); }
  private:
    using iterator_base<std_vector_base<T, strategy> >::this_;
    using iterator_base<std_vector_base<T, strategy> >::addr_;
  };

  iterator begin() {
    size_t beg = 0;
    for(; beg < std_vec_.size() && std_vec_.is_deleted(beg); ++beg);
    return iterator(beg, *this);
  };
  iterator end() { return iterator(std_vec_.size(), *this); };

  const_iterator begin() const {
    size_t beg = 0;
    for(; beg < std_vec_.size() && std_vec_.is_deleted(beg); ++beg);
    return const_iterator(beg, *this);
  }
  const_iterator end() const { return const_iterator(std_vec_.size(), *this); };

  const_iterator add(const value_type &v) {
    std_vec_.push_back(v);
    ++size_;
    return const_iterator(std_vec_.size()-1, *this);
  }

  //! @note use std::vector::erase will be very slow implementaton,
  //! and leads to non-constant iterator-element relationship.
  // void del(const const_iterator &i) { assert(0); erase(i); }

  void del(const const_iterator &i) {
    std_vec_.del(i.addr_);
    --size_;
  }
private:
  vector_t std_vec_;
  size_t size_;
};

//! @brief bit is faster in some cases.  User can choose.
template <typename T>
class std_vector : public std_vector_base<T, vector_by_bit>
{};

}}

#endif
