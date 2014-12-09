#ifndef HJ_MAT_DATA_TREE_H_
#define HJ_MAT_DATA_TREE_H_

#include <fstream>
#include <zjucad/matrix/matrix.h>

#include "../include/data_tree.h"

namespace hj { namespace data_tree {

typedef zjucad::matrix::matrix<double> matrixd_t;
typedef zjucad::matrix::matrix<size_t> matrixi_t;

template <typename TYPE>
class dt_matrix_fs_t : public io_trait
{
public:
  virtual std::string value(const boost::any &val) const {
    std::ostringstream oss;
    const matrixd_t &m = *boost::any_cast<boost::shared_ptr<zjucad::matrix::matrix<TYPE> > >(val);
    oss << m.size(1) << ", " << m.size(2);
    return oss.str();
  }
  virtual int write(boost::any &ar, const boost::any &val) const { 
    std::ofstream **ppofs = boost::any_cast<std::ofstream *>(&ar);
    if(ppofs) {
      const matrixd_t &m = *boost::any_cast<boost::shared_ptr<zjucad::matrix::matrix<TYPE> > >(val);
      ofs_archive::write_basic_type(**ppofs, m.size(1));
      ofs_archive::write_basic_type(**ppofs, m.size(2));
      (*ppofs)->write((const char *)&m[0], m.size()*sizeof(typename zjucad::matrix::matrix<TYPE>::value_type));
      return 0;
    }
    return -1;
  }
  virtual int read(boost::any &ar, boost::any &val) const {
    std::ifstream **ppifs = boost::any_cast<std::ifstream *>(&ar);
    if(ppifs) {
      size_t nrow, ncol;
      ifs_archive::read_basic_type(**ppifs, nrow);
      ifs_archive::read_basic_type(**ppifs, ncol);
      boost::shared_ptr<zjucad::matrix::matrix<TYPE> > m(new zjucad::matrix::matrix<TYPE>(nrow, ncol));
      (*ppifs)->read((char *)&(*m)[0], m->size()*sizeof(typename zjucad::matrix::matrix<TYPE>::value_type));
      val = m;
      return 0;
    }
    return -1;
  }
};

class dt_matrixd_fs_t : public dt_matrix_fs_t<matrixd_t::value_type>
{
  virtual const char* key(void) const { return "matrixd_t"; }
  virtual int write(boost::any &ar, const boost::any &val) const { 
    return dt_matrix_fs_t<matrixd_t::value_type>::write(ar, val);
  }
  virtual int read(boost::any &ar, boost::any &val) const {
    return dt_matrix_fs_t<matrixd_t::value_type>::read(ar, val);
  }
};

class dt_matrixi_fs_t : public dt_matrix_fs_t<matrixi_t::value_type>
{
  virtual const char* key(void) const { return "matrixi_t"; }
  virtual int write(boost::any &ar, const boost::any &val) const { 
    return dt_matrix_fs_t<matrixi_t::value_type>::write(ar, val);
  }
  virtual int read(boost::any &ar, boost::any &val) const {
    return dt_matrix_fs_t<matrixi_t::value_type>::read(ar, val);
  }
};

template <typename value_type>
class dt_matrix_prt_t : public io_trait
{
public:
  virtual int write(boost::any &ar, const boost::any &val) const {
    std::ostream **pos = (boost::any_cast<print_writer::ar_t *>(&ar));
    typedef typename zjucad::matrix::matrix<value_type> mat;
    const mat *ptr = boost::any_cast<const boost::shared_ptr<mat> >(val).get();
    (**pos) << '(' << ptr->size(1) << ", " << ptr->size(2) << ")\t@" << ptr;
    return 0;
  }
};

class dt_matrixd_prt_t : public dt_matrix_prt_t<double>
{
public:
  virtual const char* key(void) const { return "matrixd_t"; }
};

class dt_matrixi_prt_t : public dt_matrix_prt_t<size_t>
{
public:
  virtual const char* key(void) const { return "matrixi_t"; }
};

}}

#endif
