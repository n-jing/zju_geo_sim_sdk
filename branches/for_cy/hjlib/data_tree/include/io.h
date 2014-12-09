#ifndef HJ_DATA_TREE_IO_H_
#define HJ_DATA_TREE_IO_H_

#include "data_tree.h"

#include <fstream>

#include <boost/property_tree/ptree.hpp>

namespace hj { namespace data_tree {

class io_trait
{
public:
  virtual ~io_trait(){}
  virtual const char* key(void) const = 0;
  virtual int write(boost::any &ar, const boost::any &val) const { return -1; }
  virtual int read(boost::any &ar, boost::any &val) const { return -1; }
};

typedef std::map<const std::type_info *,
                 const std::shared_ptr<io_trait> > type_op_map;

class data_tree_reader
{
public:
  virtual ~data_tree_reader(){}
  data_tree_reader(const type_op_map &tr);
protected:
  typedef std::map<std::string, std::pair<const io_trait *, const std::type_info *> > name_to_type_t;
  name_to_type_t name_to_type_;
};

class property_tree_reader : public data_tree_reader
{
public:
  property_tree_reader(const type_op_map &tr);
  int read(const boost::property_tree::ptree &pt, data_tree &dt);
};

type_op_map default_property_tree_trait(void);

class data_tree_writer
{
public:
  virtual ~data_tree_writer(){}
  data_tree_writer(const type_op_map &tr);
protected:
  const type_op_map &tr_;
};

class print_writer : public data_tree_writer
{
public:
  typedef std::ostream ar_t;

  print_writer(const type_op_map &tr);
  int write(std::ostream &os, const data_tree &dt, size_t indent = 0);
};

type_op_map default_print_trait();

type_op_map default_fsarchive_trait(void);

class oarchive // hierarchical archive, e.g. hdf
{
public:
  virtual ~oarchive(){}
  int write(const data_tree &data, const type_op_map &type, const char *dir_name = 0);
  virtual void write_node(const char *name, const boost::any &value, const std::type_info &ti, const io_trait &tt) = 0;
};

class iarchive // hierarchical archive, e.g. hdf
{
public:
  virtual ~iarchive(){}
  iarchive(const type_op_map &type);
  int read(data_tree &data);
  virtual int read_meta(std::string &name, std::string &type_key) { return -1;}
  // NOTICE: must be boost::shared_ptr<T>, return int(-1) for not
  // supported, other for termination
  virtual int read_data(const io_trait &tt, boost::any &val) { return -1; }
private:
  std::map<std::string, const std::type_info *> name_to_type_;
  const type_op_map &type_;
};

class fs_archive
{
public:
  virtual ~fs_archive(){}
protected:
  fs_archive();
  void put_magic(std::ofstream &ofs) const;
  int status(std::ifstream &ifs) const;
private:
  const char *magic_;
};

class ofs_archive : public oarchive, public fs_archive
{
public:
  ofs_archive(const char *name);
  virtual void write_node(const char *name, const boost::any &value, const std::type_info &ti, const io_trait &tt);
  template <typename T>
  static void write_basic_type(std::ofstream &ofs, const T &val) {
    ofs.write((const char *)&val, sizeof(T));
  }
protected:
  std::ofstream ofs_;
  void write_string(const char *str);
  int write_data(const boost::any &value, const io_trait &tt);
};

class ifs_archive : public iarchive, public fs_archive
{
public:
  ifs_archive(const char *name, const type_op_map &dt);
  virtual int read_meta(std::string &name, std::string &type_key);
  virtual int read_data(const io_trait &tt, boost::any &val);
  template <typename T>
  static void read_basic_type(std::ifstream &ifs, const T &val) {
    ifs.read((char *)&val, sizeof(T));
  }
private:
  std::ifstream ifs_;
  int read_string(std::string &str);
};

}}

#endif
