#include "../include/optimizer.h"
#include "../src/math_func2hj_func.h"
#include <iostream>
#include <zjucad/matrix/io.h>
#include <hjlib/math_func/func_aux.h>
#include <hjlib/math_func/operation.h>
#include <hjlib/function/function.h>

using namespace std;
using namespace hj::math_func;

class x2_1 : public hj::math_func::math_func_t<double,int32_t>
{
public:
  virtual size_t nx() const{
    return 2;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const value_type *x, const coo2val_t<value_type, int_type> &cv,
                   func_ctx *ctx = 0) const {
    if(k == 0){
        int_type c[] = {0};
        cv[c] += *x-1;
      }
    if(k == 1){
        int_type c[] = {0,0};
        cv[c] += 1;
      }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const{
    if(k == 1){
        int_type c[] = {0,0};
        l2g.add(cs,c);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 1;
    if(k == 2) return 0;
  }
};

class sinx0: public hj::math_func::math_func_t<double,int32_t>
{
public:
  virtual size_t nx()const{
    return 1;
  }
  virtual size_t nf()const{
    return 1;
  }
  virtual int eval(size_t k, const value_type *x, const coo2val_t<value_type, int_type> &cv,
                   func_ctx *ctx = 0) const{
    if(k == 0){
        int_type c[] = {0};
        cv[c] += sin(x[0]);
      }
    if(k == 1){
        int_type c[] = {0,0};
        cv[c] += cos(x[0]);
      }
    if(k == 2){
        int_type c[] = {0,0,0};
        cv[c] += -sin(x[0]);
      }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
    if(k == 1){
        int_type c[] = {0,0};
        l2g.add(cs,c);
      }
    if(k == 2){
        int_type c[] = {0,0,0};
        l2g.add(cs,c);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 1;
    if(k == 0) return 0;
  }
};


class sinx0_cosx1: public hj::math_func::math_func_t<double,int32_t>
{
public:
  virtual size_t nx()const{
    return 2;
  }
  virtual size_t nf()const{
    return 1;
  }
  virtual int eval(size_t k, const value_type *x, const coo2val_t<value_type, int_type> &cv,
                   func_ctx *ctx = 0) const{
    if(k == 0){
        int_type c[] = {0};
        cv[c] += sin(x[0]) + cos(x[1]);
      }
    if(k == 1){
        int_type c[] = {0,0};
        cv[c] += cos(x[0]);

        int_type c2[] = {0,1};
        cv[c2] += -sin(x[1]);
      }
    if(k == 2){
        int_type c[] = {0,0,0};
        cv[c] += -sin(x[0]);
        int_type c2[] = {0,1,1};
        cv[c2] += -cos(x[1]);
      }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g, func_ctx *ctx = 0) const {
    if(k == 1){
        int_type c[] = {0,0};
        l2g.add(cs,c);

        int_type c2[] = {0,1};
        l2g.add(cs,c2);
      }
    if(k == 2){
        int_type c[] = {0,0,0};
        l2g.add(cs,c);
        int_type c2[] = {0,1,1};
        l2g.add(cs,c2);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 2;
    if(k == 2) return 2;
  }
};


int main(int argc, char * argv[])
{
  typedef math_func_t<double,int32_t> math_func_T;
  typedef shared_ptr<math_func_T> math_func_ptr;
  shared_ptr<vector<math_func_ptr> >  funcs(new vector<math_func_ptr>);
  funcs->push_back(math_func_ptr(new x2_1));
  funcs->push_back(math_func_ptr(new sinx0_cosx1));

  boost::property_tree::ptree pt;
  pt.put("package.value", "hj");
  pt.put("iter.value",10);
  zjucad::matrix::matrix<double> x = zjucad::matrix::ones<double>(2,1);

  math_func_ptr all_funcs_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(funcs));
  zjucad::optimize(*all_funcs_cat, x, pt);
  std::cerr << x << std::endl;

  return 0;
}

