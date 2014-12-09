#ifndef _TRUST_REGION_SOLVER_STEIHAUG_CG_H_
#define _TRUST_REGION_SOLVER_STEIHAUG_CG_H_

#include "trust_region_solver.h"
#include "kkt_preconditioner.h"
#include "pcg_solver.h"
#include <limits>

extern "C" {
    int solve_quadric(double c[3], double s[2]);
}

namespace zjucad
{
namespace LSCO
{

template <typename KERNEL_TYPE>
class trust_region_solver_steihaug_cg : public trust_region_solver<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_type;
    using trust_region_solver<KERNEL_TYPE>::m_opt;
    typedef trust_region_solver<KERNEL_TYPE> BASE_CLASS;
    //functions
    trust_region_solver_steihaug_cg(boost::property_tree::ptree& opt)
        :trust_region_solver<KERNEL_TYPE>(opt),
         m_sol(m_opt.add_child("trust_region_solver_pcg",boost::property_tree::ptree())) {}
    //interface
    virtual void set_hessian(boost::shared_ptr<const sparse_matrix_type> H) {
        BASE_CLASS::m_H=H;
    }
    virtual void set_hessian(boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > H) {
        BASE_CLASS::m_H_krylov=H;
    }
    virtual void set_constraint(boost::shared_ptr<const sparse_matrix_type> A) {
        BASE_CLASS::m_A=A;
    }
    virtual void assemble() {
        if(BASE_CLASS::m_H) {
            m_sol.set_A(*BASE_CLASS::m_H);
            KERNEL_TYPE::resize(KERNEL_TYPE::rows(*BASE_CLASS::m_H),m_dir);
        } else if(BASE_CLASS::m_H_krylov) {
            m_sol.set_A(BASE_CLASS::m_H_krylov);
            KERNEL_TYPE::resize(BASE_CLASS::m_H_krylov->rows(),m_dir);
        } else throw std::invalid_argument("H matrix not given");

        bool non_singular=m_opt.template get<bool>("trust_region_solver_in.nonsingular_preconditioner",false);
        boost::property_tree::ptree& pre_pt=m_opt.add_child("trust_region_solver_steihaug_cg_pre",boost::property_tree::ptree());
        m_sol.set_pre(boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >());
        if(BASE_CLASS::m_A) {
            if(m_M) {
                boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > pre(new kkt_preconditioner_M_A<KERNEL_TYPE>(pre_pt,*m_M,*BASE_CLASS::m_A,non_singular));
                m_sol.set_pre(pre);
                boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > prj(new kkt_preconditioner_schur_A<KERNEL_TYPE>(pre_pt,*BASE_CLASS::m_A,non_singular));
                BASE_CLASS::m_prj=prj;
            } else if(m_invM) {
                boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > pre(new kkt_preconditioner_invM_A<KERNEL_TYPE>(pre_pt,*m_invM,*BASE_CLASS::m_A,non_singular));
                m_sol.set_pre(pre);
                boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > prj(new kkt_preconditioner_schur_A<KERNEL_TYPE>(pre_pt,*BASE_CLASS::m_A,non_singular));
                BASE_CLASS::m_prj=prj;
            } else {
                boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > pre(new kkt_preconditioner_schur_A<KERNEL_TYPE>(pre_pt,*BASE_CLASS::m_A,non_singular));
                m_sol.set_pre(pre);
                BASE_CLASS::m_prj=pre;
            }
        } else {
            if(m_M) {
                boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > pre(new kkt_preconditioner_M<KERNEL_TYPE>(pre_pt,*m_M,non_singular));
                m_sol.set_pre(pre);
            } else if(m_invM) {
                boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > pre(new default_krylov_matrix<KERNEL_TYPE>(*m_invM));
                m_sol.set_pre(pre);
            } else if(m_invM_krylov) {
                m_sol.set_pre(m_invM_krylov);
            }
        }
        BASE_CLASS::assemble();
    }
    virtual void solve(const vector_type& g,value_type delta,vector_type& x) {
		//solve using pcg
        bool warm_start=m_opt.template get<bool>("trust_region_solver_in.warm_start",false);
        m_opt.template put<bool>("trust_region_solver_pcg.pcg_in.warm_start",warm_start);
        m_opt.template put<value_type>("trust_region_solver_pcg.pcg_in.threshold",m_opt.template get<value_type>("trust_region_solver_in.threshold",1E-5f));
        m_opt.template put<int>("trust_region_solver_pcg.pcg_in.residual_type",m_opt.template get<int>("trust_region_solver_in.residual_type",NORM_2));
        m_opt.template put<bool>("trust_region_solver_pcg.pcg_in.threshold_relative",m_opt.template get<bool>("trust_region_solver_in.threshold_relative",true));
        m_opt.template put<size_t>("trust_region_solver_pcg.pcg_in.max_iteration",m_opt.template get<size_t>("trust_region_solver_in.max_iteration",1000));
        m_opt.template put<bool>("trust_region_solver_pcg.pcg_in.trust_region",true);
        m_opt.template put<value_type>("trust_region_solver_pcg.pcg_in.trust_radius",delta);
        m_opt.template put<bool>("trust_region_solver_pcg.pcg_in.report",m_opt.template get<bool>("trust_region_solver_in.report",false));

        //safe guard
        m_opt.get_child("trust_region_solver_out").clear();
        m_opt.template put<std::string>("trust_region_solver_out.termination_type","successful");
        m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","global_minima");
        if(KERNEL_TYPE::nrm2(g) < 1E-9f) {
            KERNEL_TYPE::zero(x);
            std::cout << "already minimized" << std::endl;
            return;
        }

        //solve pcg
		m_sol.set_bound(m_bounds);
        KERNEL_TYPE::copy(g,m_dir);
        KERNEL_TYPE::scal(-1.0f,m_dir);
        m_sol.solve(m_dir,x);

        //test pcg's termination_type
        std::string pcg_output_type=m_opt.template get<std::string>("trust_region_solver_pcg.pcg_out.termination_type","max_iteration");
        if(pcg_output_type == "threshold_reach") {
            m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","global_minima");
        } else if(pcg_output_type == "negative_curvature") {
            //search along direction to boundary
            KERNEL_TYPE::copy(m_sol.get_last_direction(),m_dir);
            if(KERNEL_TYPE::dot(m_dir,m_sol.get_residual()) < 0.0f)
                KERNEL_TYPE::scal(-1.0f,m_dir);
            if(!search(delta,x,m_dir,std::numeric_limits<value_type>::max())) {
                std::cerr << "rounding error" << std::endl;
                solve_cauchy(g,delta,x);
            }
            m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","negative_curvature");
        } else if(pcg_output_type == "outside_trust_region") {
            //search along direction to boundary
            if(!search(delta,x,m_sol.get_last_direction(),1.0f)) {
                std::cerr << "rounding error" << std::endl;
                solve_cauchy(g,delta,x);
            }
            m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","on_boundary");
        } else if(pcg_output_type == "max_iteration") {
            solve_cauchy(g,delta,x);
            m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","max_iteration");
        } else throw "unknown pcg status";

        if(m_opt.template get<bool>("trust_region_solver_in.report",false))
            report(g,x);
		assert(check_bound<KERNEL_TYPE>(x,m_bounds));
	}
    //preconditioner
    virtual void set_preconditioner(boost::shared_ptr<const sparse_matrix_type> M) {
        m_M=M;
    }
    virtual void set_preconditioner_inv(boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > invM,bool sync_prj) {
        m_invM_krylov=invM;
        if(sync_prj)
            BASE_CLASS::m_prj=invM;
    }
    virtual void set_preconditioner_inv(boost::shared_ptr<const sparse_matrix_type> invM) {
        m_invM=invM;
    }
private:
    bool search(value_type delta,vector_type& x,const vector_type& dir,value_type maxV) {
        double c[3]= {
            KERNEL_TYPE::dot(x,x)-delta*delta,
            2.0f*KERNEL_TYPE::dot(x,dir),
            KERNEL_TYPE::dot(dir,dir),
        };
        double s[2];
        bool found=false;
        int nr=solve_quadric(c,s);
        for(int i=0; i<nr; i++) {
            if(s[i] >= 0.0f && s[i] <= maxV) {
				maxV=s[i];
                found=true;
                break;
            }
        }
		value_type lmt=limit_bound<KERNEL_TYPE>(x,dir,m_bounds);
		if(lmt < maxV) {
			maxV=lmt;
			found=true;
		}
		KERNEL_TYPE::axpy(maxV,dir,x);
        return found;
    }
    boost::shared_ptr<const sparse_matrix_type> m_M;
    boost::shared_ptr<const sparse_matrix_type> m_invM;
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_invM_krylov;
    //solver
    vector_type m_dir;
    pcg_solver<KERNEL_TYPE> m_sol;
};

}
}


#endif
