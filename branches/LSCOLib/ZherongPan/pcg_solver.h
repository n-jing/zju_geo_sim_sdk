#ifndef _PCG_SOLVER_H_
#define _PCG_SOLVER_H_

#include "util.h"
#include "krylov_matrix.h"
#include <boost/property_tree/ptree.hpp>

namespace zjucad
{
namespace LSCO
{

//accepted parameter,default value
//pcg_in.warm_start,false
//pcg_in.threshold,1E-5f
//pcg_in.residual_type,NORM_2
//pcg_in.threshold_relative,true
//pcg_in.max_iteration,1000
//pcg_in.trust_region,false
//pcg_in.trust_radius,none
//pcg_in.flexible,true
//pcg_in.report,false

//returned parameter
//pcg_out.iteration_count,size_t
//pcg_out.termination_type,std::string
//	threshold_reach
//	negative_curvature
//	outside_trust_region
//	max_iteration
//pcg_out.residual_[0,1,2,...]

template <typename KERNEL_TYPE>
class pcg_solver
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //functions
    pcg_solver(boost::property_tree::ptree& opt):m_opt(opt) {
        m_opt.add_child("pcg_in",boost::property_tree::ptree());
        m_opt.add_child("pcg_out",boost::property_tree::ptree());
    }
    void set_A(const sparse_matrix_type& A) {
        boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > A_krylov;
        A_krylov.reset(new default_krylov_matrix<KERNEL_TYPE>(A));
        set_A(A_krylov);
    }
    void set_A(boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > A) {
        m_A=A;
        KERNEL_TYPE::resize(m_A->rows(),m_z);
        KERNEL_TYPE::resize(m_A->rows(),m_s);
        KERNEL_TYPE::resize(m_A->rows(),m_r);
        KERNEL_TYPE::resize(m_A->rows(),m_m);
    }
    void set_pre(const sparse_matrix_type& invM) {
        boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > pre_krylov;
        pre_krylov.reset(new default_krylov_matrix<KERNEL_TYPE>(invM));
        set_pre(pre_krylov);
    }
    void set_pre(boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > pre) {
        m_pre=pre;
    }
	void set_bound(const std::vector<std::pair<size_t,value_type> >& bd){m_bounds=bd;}
    void solve(const vector_type& rhs,vector_type& result) {
        m_opt.get_child("pcg_out").clear();

        //locals
        value_type rho,rho_new,tol,alpha,beta,res;
        size_t iter,max_iter;

        //warm start option
        if(!m_opt.get<bool>("pcg_in.warm_start",false))
            KERNEL_TYPE::zero(result);

        //initialize m_r
        m_A->mul(result,m_z);
        KERNEL_TYPE::copy(rhs,m_r);
        KERNEL_TYPE::axpy(-1.0f,m_z,m_r);

        //initialize z
        if(m_pre)
            m_pre->mul(m_r,m_z);
        else KERNEL_TYPE::copy(m_r,m_z);

        //check residual
        tol=m_opt.get<value_type>("pcg_in.threshold",1E-5f);
        res=residual<KERNEL_TYPE>(m_z,m_opt.get<int>("pcg_in.residual_type",NORM_2));
        if(m_opt.get<bool>("pcg_in.threshold_relative",true))
            tol*=res;
        if(res < tol) {
            m_opt.put<size_t>("pcg_out.iteration_count",0);
            m_opt.put<std::string>("pcg_out.termination_type","threshold_reach");
            return;
        }

        //main iteration
        rho=KERNEL_TYPE::dot(m_r,m_z);
        KERNEL_TYPE::copy(m_z,m_s);

        max_iter=m_opt.get<size_t>("pcg_in.max_iteration",1000);
        for(iter=0; iter<max_iter; iter++) {
            //search along conjugate direction
            m_A->mul(m_s,m_z);
            m_curv=KERNEL_TYPE::dot(m_s,m_z);
            if(m_curv < 0.0f) {
                m_opt.put<size_t>("pcg_out.iteration_count",iter);
                m_opt.put<std::string>("pcg_out.termination_type","negative_curvature");
                return;
            }
            alpha=rho/m_curv;
			m_bound_limit=limit_bound<KERNEL_TYPE>(result,m_s,m_bounds);
			if(m_bound_limit < alpha) {
				KERNEL_TYPE::scal(alpha,m_s);
				m_opt.put<size_t>("pcg_out.iteration_count",iter);
				m_opt.put<std::string>("pcg_out.termination_type","outside_trust_region");
				return;
			}
            KERNEL_TYPE::axpy(alpha,m_s,result);
            if(m_opt.get<bool>("pcg_in.trust_region",false))
                if(KERNEL_TYPE::nrm2(result) > m_opt.get<value_type>("pcg_in.trust_radius")) {
                    //step back
                    KERNEL_TYPE::axpy(-alpha,m_s,result);
                    KERNEL_TYPE::scal(alpha,m_s);
                    m_opt.put<size_t>("pcg_out.iteration_count",iter);
                    m_opt.put<std::string>("pcg_out.termination_type","outside_trust_region");
                    return;
                }
            if(m_opt.get<bool>("pcg_in.flexible",true))
                KERNEL_TYPE::copy(m_r,m_m);

            //update residual and check
            KERNEL_TYPE::axpy(-alpha,m_z,m_r);
            if(m_pre)
                m_pre->mul(m_r,m_z);
            else KERNEL_TYPE::copy(m_r,m_z);
            res=residual<KERNEL_TYPE>(m_z,m_opt.get<int>("pcg_in.residual_type",NORM_2));
            if(res < tol) {
                m_opt.put<size_t>("pcg_out.iteration_count",iter);
                m_opt.put<std::string>("pcg_out.termination_type","threshold_reach");
                return;
            }
            if(m_opt.get<bool>("pcg_in.report",false)) {
                std::ostringstream oss;
                oss << "pcg_out.residual_" << iter;
                m_opt.put<value_type>(oss.str(),res);
            }

            //generate new search direction
            rho_new=KERNEL_TYPE::dot(m_z,m_r);
            if(m_opt.get<bool>("pcg_in.flexible",true))
                beta=(rho_new-KERNEL_TYPE::dot(m_z,m_m))/rho;
            else beta=rho_new/rho;
            KERNEL_TYPE::axpy(beta,m_s,m_z);
            KERNEL_TYPE::copy(m_z,m_s);
            rho=rho_new;
        }

        m_opt.put<size_t>("pcg_out.iteration_count",iter);
        m_opt.put<std::string>("pcg_out.termination_type","max_iteration");
        return;
    }
    value_type get_last_curv() const {
        return m_curv;
    }
	value_type get_bound_limit() const{return m_bound_limit;}
    const vector_type& get_last_direction() const {
        return m_s;
    }
    const vector_type& get_residual() const {
        return m_r;
    }
    void report() {
        size_t nr_iter=m_opt.get<size_t>("pcg_out.iteration_count",0);
        std::cout << "iter count: " << nr_iter << std::endl;
        std::cout << "termination type: " << m_opt.get<std::string>("pcg_out.termination_type","threshold_reach") << std::endl;
        if(m_opt.get<bool>("pcg_in.report",false)) {
            for(size_t iter=0; iter<nr_iter; iter++) {
                std::ostringstream oss;
                oss << "pcg_out.residual_" << iter;
                std::cout << "\tResidual " << iter << ": " << m_opt.get<value_type>(oss.str(),0.0f) << std::endl;
            }
        }
        std::cout << std::endl;
    }
private:
    value_type m_curv,m_bound_limit;
    vector_type m_z,m_s,m_r,m_m;
	std::vector<std::pair<size_t,value_type> > m_bounds;
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_A;
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_pre;
    boost::property_tree::ptree& m_opt;
};

}
}

#endif
