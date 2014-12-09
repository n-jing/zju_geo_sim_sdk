#ifndef _TRUST_REGION_SOLVER_DOUBLE_DOGLEG_SOLVER_H_
#define _TRUST_REGION_SOLVER_DOUBLE_DOGLEG_SOLVER_H_

#include "trust_region_solver.h"
#include "kkt_preconditioner.h"

extern "C" {
    int solve_quadric(double c[3], double s[2]);
}

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename MATRIX_TYPE>
struct direct_solver_kernel;

///
/// @brief trust region solver with double dogleg search, needs gamma as input
///
template <typename KERNEL_TYPE>
class trust_region_solver_double_dogleg : public trust_region_solver<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_kernel_type;
    typedef trust_region_solver<KERNEL_TYPE> BASE_CLASS;
    //functions
    trust_region_solver_double_dogleg(boost::property_tree::ptree& opt)
        :trust_region_solver<KERNEL_TYPE>(opt) {}
	virtual void set_hessian(boost::shared_ptr<const sparse_matrix_type> H) {
        BASE_CLASS::m_H=H;
    }
    virtual void set_hessian_ATA(const sparse_matrix_type& A) {
        BASE_CLASS::m_H_krylov.reset(new krylov_matrix_ATA<KERNEL_TYPE>(A));
        bool non_singular=BASE_CLASS::m_opt.template get<bool>("trust_region_solver_in.nonsingular_preconditioner",false);
        boost::property_tree::ptree& pre_pt=BASE_CLASS::m_opt.add_child("trust_region_solver_double_dogleg_pre",boost::property_tree::ptree());
        m_sol.reset(new kkt_preconditioner_ATA<KERNEL_TYPE>(pre_pt,A,non_singular));
    }
    virtual void set_constraint(boost::shared_ptr<const sparse_matrix_type> A) {
        BASE_CLASS::m_A = A;
    }
    virtual void assemble() {
        bool non_singular=BASE_CLASS::m_opt.template get<bool>("trust_region_solver_in.nonsingular_preconditioner",false);
        boost::property_tree::ptree& pre_pt=BASE_CLASS::m_opt.add_child("trust_region_solver_double_dogleg_pre",boost::property_tree::ptree());
        if(BASE_CLASS::m_H_krylov) {
            if(BASE_CLASS::m_H)
                throw "explicit hessian not supported in ATA or AAT mode";
            if(BASE_CLASS::m_A)
                throw "constraint not supported in ATA or AAT mode";
            BASE_CLASS::m_prj.reset((krylov_matrix<KERNEL_TYPE>*)NULL);
        } else if(BASE_CLASS::m_A) {
            m_sol.reset(new kkt_preconditioner_M_A<KERNEL_TYPE>(pre_pt,*BASE_CLASS::m_H,*BASE_CLASS::m_A,non_singular));
            BASE_CLASS::m_prj.reset(new kkt_preconditioner_schur_A<KERNEL_TYPE>(pre_pt,*BASE_CLASS::m_A,non_singular));
        } else {
            m_sol.reset(new kkt_preconditioner_M<KERNEL_TYPE>(pre_pt,*BASE_CLASS::m_H,true));
            BASE_CLASS::m_prj.reset((krylov_matrix<KERNEL_TYPE>*)NULL);
        }

        if(BASE_CLASS::m_H) {
            KERNEL_TYPE::resize(KERNEL_TYPE::rows(*BASE_CLASS::m_H),m_dir);
            KERNEL_TYPE::resize(KERNEL_TYPE::rows(*BASE_CLASS::m_H),m_gc);
            KERNEL_TYPE::resize(KERNEL_TYPE::rows(*BASE_CLASS::m_H),m_gn);
        } else if(BASE_CLASS::m_H_krylov) {
            KERNEL_TYPE::resize(BASE_CLASS::m_H_krylov->rows(),m_dir);
            KERNEL_TYPE::resize(BASE_CLASS::m_H_krylov->rows(),m_gc);
            KERNEL_TYPE::resize(BASE_CLASS::m_H_krylov->rows(),m_gn);
        } else throw std::invalid_argument("H matrix not given");
        BASE_CLASS::assemble();
    }
    virtual void solve(const vector_type& g,value_type delta,vector_type& x) {
        //safe guard
        BASE_CLASS::m_opt.get_child("trust_region_solver_out").clear();
        BASE_CLASS::m_opt.template put<std::string>("trust_region_solver_out.termination_type","successful");
        BASE_CLASS::m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","global_minima");
        if(KERNEL_TYPE::nrm2(g) < 1E-9f) {
            KERNEL_TYPE::zero(x);
            //std::cout << "already minimized" << std::endl;
            return;
        }

        //first phase
		value_type min_f;
        bool block=false;
        {
            block=solve_cauchy(g,delta,m_gc);
            KERNEL_TYPE::copy(m_gc,x);
            min_f=eval(g,m_gc);
            if(block)
                BASE_CLASS::m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","on_boundary_seg0");
        }

        //calculate newton point
        KERNEL_TYPE::copy(g,m_dir);
        KERNEL_TYPE::scal(-1.0f,m_dir);
        m_sol->mul(m_dir,m_gn);

        //second phase
        value_type denom=KERNEL_TYPE::dot(m_gc,m_gn);
        if(!block && std::abs(denom) > 1E-6f) {
            //decide optimal gamma
            value_type gamma_star=KERNEL_TYPE::dot(m_gc,m_gc)/denom;
            KERNEL_TYPE::copy(m_gn,m_dir);
            KERNEL_TYPE::scal(gamma_star,m_dir);
            assert(eval(g,m_dir) <= eval(g,m_gc)+1E-6f);
            assert(KERNEL_TYPE::nrm2(m_dir)+1E-6f >= KERNEL_TYPE::nrm2(m_gc));
            KERNEL_TYPE::axpy(-1.0f,m_gc,m_dir);

            // find range which statisfiy |x| < delta
            value_type max_tau=std::min<value_type>(1.0f,limit_bound<KERNEL_TYPE>(m_gc,m_dir,BASE_CLASS::m_bounds));
            {
                value_type c[3]= {
                    KERNEL_TYPE::dot(m_gc,m_gc)-delta*delta,
                    2*KERNEL_TYPE::dot(m_gc,m_dir),
                    KERNEL_TYPE::dot(m_dir,m_dir)
                };
                value_type roots[2];
                int nr=solve_quadric(c,roots);
                for(int i=0; i<nr; i++)
                    if(roots[i] > 0.0f && roots[i] < max_tau) {
                        max_tau=roots[i];
                        block=true;
                    }
            }
            KERNEL_TYPE::axpy(max_tau,m_dir,m_gc);
            value_type f=eval(g,m_gc);
            if(f < min_f) {
                min_f=f;
                KERNEL_TYPE::copy(m_gc,x);
            }
            if(block)
                BASE_CLASS::m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","on_boundary_seg1");
        }
        //third phase
        if(!block) {
            value_type norm_n=KERNEL_TYPE::nrm2(m_gn);
            value_type thres=delta/norm_n;
                        value_type thres_lmt=limit_bound<KERNEL_TYPE>(x,m_gn,BASE_CLASS::m_bounds);
            if(norm_n < 1E-6f || thres > 1.0f && thres_lmt > 1.0f) {
                thres=1.0f;
                BASE_CLASS::m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","global_minima");
            } else BASE_CLASS::m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","on_boundary_seg2");
            KERNEL_TYPE::copy(m_gn,x);
            KERNEL_TYPE::scal(std::min(thres,thres_lmt),x);
        }

        if(BASE_CLASS::m_opt.template get<bool>("trust_region_solver_in.report",false))
            report(g,x);
                assert(check_bound<KERNEL_TYPE>(x,BASE_CLASS::m_bounds));
    }
protected:
    //solver
    vector_type m_dir,m_gn,m_gc;
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_sol;
};

}
}

#endif
