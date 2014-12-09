#ifndef HJ_DATA_TREE_H_
#define HJ_DATA_TREE_H_

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <stdexcept>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

namespace hj { namespace data_tree {

class data_tree
{
public:
  typedef std::map<std::string, std::pair<boost::any, const std::type_info *> > container;
  typedef container::const_iterator const_iterator;
  
  const_iterator begin(void) const { return dat_.begin(); }
  const_iterator end(void) const { return dat_.end(); }

  //! NOTICE: ptr must be a heap ptr, the ownership transfers to dt.
  template<typename T>
  boost::shared_ptr<T> put(const char *path, T *ptr) {
    std::pair<data_tree *, const char *> db = get_or_create_dir_base_name(path);
    boost::shared_ptr<T> rtn(ptr);
    db.first->dat_[db.second] = make_pair(rtn, &typeid(T));
    return rtn;
  }

  template<typename T>
  boost::shared_ptr<T> put(const char *path, const T &v) {
    return put(path, new T(v));
  }

  ///! can be used to create link
  template<typename T>
  boost::shared_ptr<T> put(const char *path, const boost::shared_ptr<T> &v) {
    std::pair<data_tree *, const char *> db = get_or_create_dir_base_name(path);
    db.first->dat_[db.second] = make_pair(v, &typeid(T));
    return v;
  }

  template<typename T>
  boost::shared_ptr<T> get(const char *path, const char *FILE = "", size_t LINE = 0) {
    boost::any *entry = get_entry(path);
    report(entry, path, FILE, LINE);
    return boost::any_cast<boost::shared_ptr<T> >(*entry);
  }

  template<typename T>
  const boost::shared_ptr<const T> get(const char *path, const char *FILE = 0, size_t LINE = 0) const {
    const boost::any *entry = get_entry(path);
    report(entry, path, FILE, LINE);
    return boost::any_cast<const boost::shared_ptr<T> >(*entry);
  }

  template <typename T>
  const T get_opt_val(const char *path, const T &opt, const char *FILE = 0, size_t LINE = 0) const {
    const boost::any *entry = get_entry(path);
    if(entry) return *boost::any_cast<const boost::shared_ptr<T> >(*entry);
    return opt;
  }

  bool has(const char *path) const {
    return get_entry(path) != 0;
  }

  const container &get(void) const { return dat_; }

  // ptr must be boost::shared_ptr<T>, ti is typeid(T)
  void put(const char *path, const boost::any &ptr, const std::type_info *ti);
private:
  boost::any *get_entry(const char *path);
  const boost::any *get_entry(const char *path) const;
  void report(const boost::any *entry, const char *path, const char *FILE, size_t LINE) const;
  /// @brief: tree[path] = tree[dirname][base] = rtn.first[rtn.second], assuming
  /// simple backslash format relative path
  std::pair<const data_tree *, const char *> get_dir_base_name(const char *path) const;
  std::pair<data_tree *, const char *> get_dir_base_name(const char *path);
  std::pair<data_tree *, const char *> get_or_create_dir_base_name(const char *path);
  container dat_;
};

}}

#endif
