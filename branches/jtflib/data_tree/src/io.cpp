#include "../include/io.h"

#include <iostream>
#include <stdexcept>

#include <zjucad/matrix/matrix.h>

using namespace std;

namespace hj { namespace data_tree {

typedef zjucad::matrix::matrix<size_t> matrixi_t;
typedef zjucad::matrix::matrix<double> matrixd_t;

#define DEF_FSARCHIVE_TRAIT(TYPE)                                       \
  class dt_fsarchive_##TYPE : public io_trait                           \
  {                                                                     \
  public:                                                               \
    virtual const char* key(void) const { return #TYPE; }               \
    virtual int write(boost::any &ar, const boost::any &val) const {    \
      ofstream **ppofs = boost::any_cast<ofstream *>(&ar);              \
      if(ppofs) {                                                       \
        const TYPE &v = *boost::any_cast<boost::shared_ptr<TYPE> >(val); \
        (*ppofs)->write((const char *)&v, sizeof(TYPE));                \
        return 0;                                                       \
      }                                                                 \
      return -1;                                                        \
    }                                                                   \
    virtual int read(boost::any &ar, boost::any &val) const {           \
      ifstream **ppifs = boost::any_cast<ifstream *>(&ar);              \
      if(ppifs) {                                                       \
        TYPE v;                                                         \
        (*ppifs)->read((char *)&v, sizeof(TYPE));                       \
        val = boost::shared_ptr<TYPE>(new TYPE(v));                     \
        return 0;                                                       \
      }                                                                 \
      return -1;                                                        \
    }                                                                   \
  };

#define INS_FSARCHIVE_TRAIT(TYPE)                                 \
    rtn.insert(make_pair(&typeid(TYPE),  std::shared_ptr<hj::data_tree::io_trait>(new dt_fsarchive_##TYPE)));

#define ADD_FSARCHIVE_TRAIT(TYPE)               \
  DEF_FSARCHIVE_TRAIT(TYPE);                    \
  INS_FSARCHIVE_TRAIT(TYPE);

class dt_fsarchive_data_tree : public io_trait
{
  virtual const char* key(void) const { return "data_tree"; }
};

#define ADD_BASIC_TYPES(MACRO)                  \
    MACRO(bool);                                \
    MACRO(char);                                \
    MACRO(int);                                 \
    MACRO(long);                                \
    MACRO(float);                               \
    MACRO(double);                              \
    MACRO(size_t);                              \

type_op_map default_fsarchive_trait(void)
{
  type_op_map rtn;
  ADD_BASIC_TYPES(ADD_FSARCHIVE_TRAIT);
  INS_FSARCHIVE_TRAIT(data_tree);
  return rtn;
}

int oarchive::write(const data_tree &data, const type_op_map &type, const char *dir_name)
{
  std::string cur_dir;
  if(dir_name) cur_dir = dir_name;

  for(data_tree::const_iterator di = data.begin(); di != data.end(); ++di) {
    const type_op_map::const_iterator ti = type.find(di->second.second);
    if(ti == type.end()) continue; // no trait to serialize
    io_trait &tt = *ti->second;
    std::string path = cur_dir+(cur_dir.empty()?"":"/")+di->first;
    if(di->second.second == &typeid(data_tree)) { // group
      const data_tree &sub_data = *data.get<data_tree>(di->first.c_str());
      write(sub_data, type, path.c_str());
    }
    else {
      write_node(path.c_str(), di->second.first, *di->second.second, tt); // value
    }
  }
  return 0;
}

iarchive::iarchive(const type_op_map &type)
  :type_(type) {
  for(type_op_map::const_iterator i = type_.begin(); i != type_.end(); ++i)
    name_to_type_.insert(make_pair(i->second->key(), i->first));
}

int iarchive::read(data_tree &data)
{
  while(1) {
    std::string name, type_key;
    if(read_meta(name, type_key))
      break;
    const std::map<std::string, const std::type_info *>::const_iterator ti
      = name_to_type_.find(type_key);
    if(ti == name_to_type_.end()) {
      std::cout << type_key << std::endl;
      for(std::map<std::string, const std::type_info *>::const_iterator
            i = name_to_type_.begin(); i != name_to_type_.end(); ++i) {
        std::cout << i->first << std::endl;
      }
      throw std::runtime_error("unknown data type");
    }

    if(ti->second == &typeid(data_tree)) {
      data_tree &sub_data = *data.put(name.c_str(), new data_tree());
      read(sub_data);
    }
    else {
      type_op_map::const_iterator tti = type_.find(ti->second);
      assert(tti != type_.end());
      boost::any val;
      const int err = read_data(*tti->second, val);
      if(!err)
        data.put(name.c_str(), val, ti->second);
      else  {
        if(err == -1)
          continue;
        return 1;
      }
    }
  }
  return 0;
}

fs_archive::fs_archive():magic_("node_beg") {}

void fs_archive::put_magic(ofstream &ofs) const
{
  ofs.write(magic_, strlen(magic_));
}

int fs_archive::status(ifstream &ifs) const
{
  std::string magic(magic_);
  ifs.read(&magic[0], strlen(magic_));
  if(!ifs.good()) return 1;
  if(magic_ != magic_) return -1;
  return 0;
}

ofs_archive::ofs_archive(const char *name)
  :ofs_(name, std::ofstream::binary) {
}

void ofs_archive::write_node(const char *name, const boost::any &value, const std::type_info &ti, const io_trait &tt) {
  put_magic(ofs_);
  write_string(name);
  write_string(tt.key());
  write_data(value, tt);
}

void ofs_archive::write_string(const char *str)
{
  size_t len = strlen(str);
  write_basic_type(ofs_, len);
  ofs_.write(str, len);
}

int ofs_archive::write_data(const boost::any &value, const io_trait &tt)
{
  boost::any ar(&ofs_);
  return tt.write(ar, value);
}

ifs_archive::ifs_archive(const char *name, const type_op_map &dt)
  :iarchive(dt), ifs_(name, std::ifstream::binary) {
}

int ifs_archive::read_meta(std::string &name, std::string &type_key)
{
  int st = status(ifs_);
  if(st) return st;
  read_string(name);
  read_string(type_key);
  return 0;
}

// TODO: must be boost::shared_ptr<T>
int ifs_archive::read_data(const io_trait &tt, boost::any &val)
{
  boost::any ar(&ifs_);
  return tt.read(ar, val);
}

int ifs_archive::read_string(std::string &str)
{
  size_t len;
  read_basic_type(ifs_, len);
  if(ifs_.fail())
    return 1;
  str.resize(len);
  ifs_.read(&str[0], len);
  return 0;
}

data_tree_reader::data_tree_reader(const type_op_map &tr)
{
  for(type_op_map::const_iterator i = tr.begin(); i != tr.end(); ++i)
    name_to_type_.insert(make_pair(i->second->key(), make_pair(i->second.get(), i->first)));
}

property_tree_reader::property_tree_reader(const type_op_map &tr)
  :data_tree_reader(tr)
{
}

int property_tree_reader::read(const boost::property_tree::ptree &pt, data_tree &dt)
{
  boost::property_tree::ptree::const_iterator i;
  for(i = pt.begin(); i != pt.end(); ++i) {
    if(i->second.begin() != i->second.end())
      read(i->second, *dt.put(i->first.c_str(), new data_tree()));
    else {
      string word;
      istringstream iss(i->first.c_str());
      iss >> word;
      name_to_type_t::const_iterator t
        = name_to_type_.find(word);
      if(t == name_to_type_.end()) {
        cerr << "# unkown type: " << word << endl;
        continue;
      }
      boost::any val;
      boost::any ar(i->second.get_value<string>().c_str());
      t->second.first->read(ar, val);
      iss >> word;
      dt.put(word.c_str(), val, t->second.second);
    }
  }
  return 0;
}

#define DEF_PROPERTY_TREE_TRAIT(TYPE)                         \
  class property_tree_##TYPE : public io_trait                \
  {                                                           \
  public:                                                     \
    virtual const char* key(void) const { return #TYPE; }     \
    virtual int read(boost::any &ar, boost::any &val) const { \
      const char **str = boost::any_cast<const char *>(&ar);  \
      TYPE r;                                                 \
      istringstream iss(*str);                                \
      iss >> r;                                               \
      val = boost::shared_ptr<TYPE>(new TYPE(r));             \
      return 0;                                               \
    }                                                         \
  };

#define INS_PROPERTY_TREE_TRAIT(TYPE)                             \
    rtn.insert(make_pair(&typeid(TYPE),  std::shared_ptr<hj::data_tree::io_trait>(new property_tree_##TYPE)));

#define ADD_PROPERTY_TREE_TRAIT(TYPE)           \
  DEF_PROPERTY_TREE_TRAIT(TYPE);                \
  INS_PROPERTY_TREE_TRAIT(TYPE);

class property_tree_string : public io_trait
{
public:
  virtual const char* key(void) const { return "string"; }
  virtual int read(boost::any &ar, boost::any &val) const {
    const char **str = boost::any_cast<const char *>(&ar);
    string r(*str);
    val = boost::shared_ptr<string>(new string(r));
    return 0;
  }
};

type_op_map default_property_tree_trait(void)
{
  type_op_map rtn;
  ADD_BASIC_TYPES(ADD_PROPERTY_TREE_TRAIT);
  INS_PROPERTY_TREE_TRAIT(string);
  return rtn;
}

data_tree_writer::data_tree_writer(const type_op_map &tr)
  :tr_(tr)
{}

print_writer::print_writer(const type_op_map &tr)
  :data_tree_writer(tr)
{}

int print_writer::write(ostream &os, const data_tree &dt, size_t indent)
{
  boost::any ar(&os);
  for(data_tree::const_iterator di = dt.begin(); di != dt.end(); ++di) {
    for(size_t i = 0; i < indent; ++i) os << "    ";
    const type_op_map::const_iterator ti = tr_.find(di->second.second);
    if(ti != tr_.end()) {
      const io_trait &tt = *ti->second;
      os << di->first << ":\t" << tt.key() << '\t';
      tt.write(ar, di->second.first);
      os << '\n';
    }
    else
      os << di->first << ": " << di->second.second->name() << "\n";
    if(di->second.second == &typeid(data_tree)) {
      const data_tree &sub_dt = *dt.get<data_tree>(di->first.c_str());
      write(os, sub_dt, indent+1);
    }
  }
  return 0;
}

#define DEF_PRINT_TRAIT(TYPE)                         \
  class print_##TYPE : public io_trait                \
  {                                                           \
  public:                                                     \
    virtual const char* key(void) const { return #TYPE; }     \
    virtual int write(boost::any &ar, const boost::any &val) const { \
      ostream **pos = (boost::any_cast<print_writer::ar_t *>(&ar));  \
      const TYPE *ptr = boost::any_cast<const boost::shared_ptr<TYPE> >(val).get(); \
      (**pos) << *ptr << "\t@" << ptr;                                  \
      return 0;                                               \
    }                                                         \
  };

#define INS_PRINT_TRAIT(TYPE)                             \
    rtn.insert(make_pair(&typeid(TYPE),  std::shared_ptr<hj::data_tree::io_trait>(new print_##TYPE)));

#define ADD_PRINT_TRAIT(TYPE)           \
  DEF_PRINT_TRAIT(TYPE);                \
  INS_PRINT_TRAIT(TYPE);

ostream &operator<<(ostream &os, const data_tree dt) {
  return os;
}

type_op_map default_print_trait(void)
{
  type_op_map rtn;
  ADD_BASIC_TYPES(ADD_PRINT_TRAIT);
  ADD_PRINT_TRAIT(string);
  ADD_PRINT_TRAIT(data_tree);
  return rtn;
}

}}
