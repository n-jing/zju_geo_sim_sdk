#ifndef ZJUCAD_OPTIMIZER_JTF_H_
#define ZJUCAD_OPTIMIZER_JTF_H_

#include <boost/property_tree/ptree.hpp>
#include <jtflib/function/function.h>
#include <zjucad/matrix/matrix.h>
#include <memory>
#include "opt.h"
namespace jtf {

  int optimize(
      jtf::function::functionN1_t<double,int32_t> &f,
      zjucad::matrix::matrix<double> &x,
      boost::property_tree::ptree & pt,
      const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > *eqn_cons,
      const std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > *ineqn_cons,
      jtf::opt::callbacks * cb);

}

#endif
