#ifndef _DL_PTESC_H_
#define _DL_PTESC_H_

#include <boost/property_tree/ptree.hpp>

class PETsc
{
public:
  virtual ~PETsc() {}
};// pets_init;

class PETsc_CG
{
public:
  virtual int solve(const double *b, double *x, size_t rhs, boost::property_tree::ptree &opts) = 0;
  virtual ~PETsc_CG() {}
};

#endif /*_DL_PTESC_H_*/

