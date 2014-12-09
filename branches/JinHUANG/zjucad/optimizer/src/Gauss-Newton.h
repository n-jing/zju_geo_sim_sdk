#ifndef HJ_GAUSS_NEWTON_H_
#define HJ_GAUSS_NEWTON_H_

#include <stdint.h>

//! some of the function are not copyable, so pointer is a good choice
//! for these functions.
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>

#include <vector>

//#include <hjlib/arg_opts/arg_opts.h>
#include <hjlib/function/function.h>
#include <zjucad/linear_solver/linear_solver.h>

class nlsq_optimizer
{
public:
	virtual ~nlsq_optimizer(){}
	virtual int iterate(double *x,  const hj::function::function &f, double *residual, boost::property_tree::ptree &opts) = 0;
};

class Gauss_Newton : public nlsq_optimizer
{
public:
	Gauss_Newton(const hj::function::function &f); // used to query working space size
	virtual int iterate(double *x, const hj::function::function &f, double *residual, boost::property_tree::ptree &opts);
protected:
	int compute_g(const double *x, const hj::function::function &f, double *residual, zjucad::matrix::matrix<double> &g,
				  double *norm2 = 0);// norm is a length 2 vector
	int solve_step(const double *x, double *s, boost::property_tree::ptree &opts);
	hj::sparse::csc<double, int32_t> JT_, JTJ_;
	zjucad::matrix::matrix<double> g_, s_;
	std::auto_ptr<linear_solver> slv_;
	bool is_JT_sorted_;

	std::auto_ptr<hj::function::func_ctx> ctx_;
};

class damped_Gauss_Newton : public Gauss_Newton
{
public:
	damped_Gauss_Newton(const hj::function::function &f); // used to query working space size
	virtual int iterate(double *x, const hj::function::function &f, double *residual, boost::property_tree::ptree &opts);
};

class LM_Gauss_Newton : public Gauss_Newton
{
public:
	LM_Gauss_Newton(const hj::function::function &f);
	virtual int iterate(double *x, const hj::function::function &f, double *residual, boost::property_tree::ptree& pt);
protected:
	virtual int solve_step(const double *x, double mu, double *s, boost::property_tree::ptree &opts);
	virtual double adjust_mu(double mu, double ratio);
	virtual void compute_D(double *D);
};

class More_Gauss_Newton : public LM_Gauss_Newton
{
public:
	More_Gauss_Newton(const hj::function::function &f);
protected:
	virtual int solve_step(const double *x, double mu, double *s, boost::property_tree::ptree &opts);
//	virtual double adjust_mu(double mu, double ratio); // Seems too complex
	virtual void compute_D(double *D);
	zjucad::matrix::matrix<double> D_;
};

#endif
