#ifndef _TRUST_REGION_SOLVER_H_
#define _TRUST_REGION_SOLVER_H_

#include <limits>
#include "util.h"
#include "krylov_matrix.h"
#include "direct_solver_kernel_report.h"
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>

namespace zjucad
{
namespace LSCO
{

//accepted parameter,default value
//trust_region_solver_in.warm_start,false						//only used by trust_region_solver_steihuag_cg
//trust_region_solver_in.threshold,1E-5f						//only used by trust_region_solver_steihuag_cg
//trust_region_solver_in.residual_type,NORM_2					//only used by trust_region_solver_steihuag_cg
//trust_region_solver_in.threshold_relative,true				//only used by trust_region_solver_steihuag_cg
//trust_region_solver_in.max_iteration,1000						//only used by trust_region_solver_steihuag_cg
//trust_region_solver_in.nonsingular_preconditioner,false		//only used by trust_region_solver_steihuag_cg
//trust_region_solver_in.report,false

//returned parameter
//trust_region_solver_out.termination_type_detailed,std::string
//	global_minima
//	negative_curvature
//	on_boundary
//	on_boundary_seg0
//	on_boundary_seg1
//	on_boundary_seg2
//  max_iteration
//trust_region_solver_out.termination_type,std::string
//	successful
//	fail

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename KERNEL_TYPE>
class trust_region_solver
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //functions
    trust_region_solver(boost::property_tree::ptree& opt):m_opt(opt) {
        m_opt.add_child("trust_region_solver_in",boost::property_tree::ptree());
        m_opt.add_child("trust_region_solver_out",boost::property_tree::ptree());
    }
    virtual ~trust_region_solver() {}
    virtual void set_hessian(boost::shared_ptr<const sparse_matrix_type> H) {
        throw "not supported!";
    }
    virtual void set_hessian(boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > H) {
        throw "not supported!";
    }
    virtual void set_constraint(boost::shared_ptr<const sparse_matrix_type> A) {
        throw "not supported!";
    }
    virtual void set_constraint(boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > A) {
        throw "not supported!";
    }
	virtual void clear_bound(){m_bounds.clear();}
	virtual void set_bound(size_t off,value_type min_value){
		m_bounds.push_back(std::pair<size_t,value_type>(off,min_value));
	}
    virtual void assemble() {
        if(m_H) {
            KERNEL_TYPE::resize(KERNEL_TYPE::rows(*m_H),m_tmp);
        } else if(m_H_krylov) {
            KERNEL_TYPE::resize(m_H_krylov->rows(),m_tmp);
        } else throw std::invalid_argument("H matrix not given");
    }
    virtual void solve(const vector_type& g,value_type delta,vector_type& x)=0;
    void report(const vector_type& g,const vector_type& x) const {
        std::cout << "termination_type: " << m_opt.template get<std::string>("trust_region_solver_out.termination_type") <<  std::endl;
        std::cout << "detailed: " << m_opt.template get<std::string>("trust_region_solver_out.termination_type_detailed") << std::endl;
    }
protected:
    value_type eval(const vector_type& g,const vector_type& x) {
        if(m_H) {
            KERNEL_TYPE::mul(*m_H,x,m_tmp);
        } else if(m_H_krylov) {
            m_H_krylov->mul(x,m_tmp);
        } else throw std::invalid_argument("H matrix not given");
        return KERNEL_TYPE::dot(m_tmp,x)*0.5f+KERNEL_TYPE::dot(g,x);
    }
    //return if it's on boundary or blocked by user bound
    bool solve_cauchy(const vector_type& g,value_type delta,vector_type& x) {
        //force g to satisfy constraint
        if(m_prj)
            m_prj->mul(g,x);
        else KERNEL_TYPE::copy(g,x);
		value_type max_c=limit_bound_n<KERNEL_TYPE>(x,m_bounds);

        //search
        value_type ng=KERNEL_TYPE::nrm2(x);
        value_type curv;
        if(m_H) {
            KERNEL_TYPE::mul(*m_H,x,m_tmp);
        } else if(m_H_krylov) {
            m_H_krylov->mul(x,m_tmp);
        } else throw std::invalid_argument("H matrix not given");
        curv=KERNEL_TYPE::dot(x,m_tmp);
        if(curv <= 0.0f || ng*ng*ng > delta*curv) {
            KERNEL_TYPE::scal(-std::min(delta/ng,max_c),x);
            return true;
        } else {
            KERNEL_TYPE::scal(-std::min(ng*ng/curv,max_c),x);
            return ng*ng/curv > max_c;
        }
    }
    //input
    boost::property_tree::ptree& m_opt;
    boost::shared_ptr<const sparse_matrix_type> m_H;
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_H_krylov;
    boost::shared_ptr<const sparse_matrix_type> m_A;
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_prj;
	std::vector<std::pair<size_t,value_type> > m_bounds;
private:
    vector_type m_tmp;
};

}
}

#endif