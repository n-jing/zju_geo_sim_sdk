#include "../include/data_tree.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace hj { namespace data_tree {

void data_tree::put(const char *path, const boost::any &ptr, const std::type_info *ti)
{
    std::pair<data_tree *, const char *> db = get_or_create_dir_base_name(path);
    db.first->dat_[db.second] = make_pair(ptr, ti);
}

boost::any *data_tree::get_entry(const char *path)
{
  std::pair<data_tree *, const char *> db = get_dir_base_name(path);
  if(!db.first) return 0;
  container::iterator pos = db.first->dat_.find(db.second);
  if(pos == db.first->dat_.end()) return 0;
  return &pos->second.first;
}

void data_tree::report(const boost::any *entry, const char *path, const char *FILE, size_t LINE) const
{
  if(!entry) {
    std::ostringstream oss;
    oss << FILE << ':' << LINE << ": data_tree::get(" << path << ')';
    std::cout << oss.str() << std::endl;
    throw std::runtime_error(oss.str().c_str());
  }
}

const boost::any *data_tree::get_entry(const char *path) const
{
  std::pair<const data_tree *, const char *> db = get_dir_base_name(path);
  if(!db.first) return 0;
  container::const_iterator pos = db.first->dat_.find(db.second);
  if(pos == db.first->dat_.end()) return 0;
  return &pos->second.first;
}
  
std::pair<const data_tree *, const char *> data_tree::get_dir_base_name(const char *path) const {
  const char *beg, *end;
  std::string name;
  const data_tree *cur = this;
  for(beg = path, end = strchr(beg, '/'); cur && end; beg = end+1, end = strchr(beg, '/')) {
    name.assign(beg, end-beg);
    cur = cur->get<data_tree>(name.c_str()).get();
  }
  return std::make_pair(cur, beg);
}

std::pair<data_tree *, const char *> data_tree::get_dir_base_name(const char *path)
{
  const char *beg, *end;
  std::string name;
  data_tree *cur = this;
  for(beg = path, end = strchr(beg, '/'); cur && end; beg = end+1, end = strchr(beg, '/')) {
    name.assign(beg, end-beg);
    cur = cur->get<data_tree>(name.c_str()).get();
  }
  return std::make_pair(cur, beg);
}

std::pair<data_tree *, const char *> data_tree::get_or_create_dir_base_name(const char *path)
{
  const char *beg, *end;
  std::string name;
  data_tree *cur = this;
  for(beg = path, end = strchr(beg, '/'); cur && end; beg = end+1, end = strchr(beg, '/')) {
    name.assign(beg, end-beg);
    if(cur->has(name.c_str())) {
      cur = cur->get<data_tree>(name.c_str()).get();
    }
    else
      cur = cur->put(name.c_str(), new data_tree()).get();
  }
  return std::make_pair(cur, beg);
}

}}
