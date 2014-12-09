#ifndef _LSCO_SOLVER_INTERIOR_CG_H_
#define _LSCO_SOLVER_INTERIOR_CG_H_

#include <limits>
#include "util.h"
#include "objective_function.h"
#include "trust_region_solver_double_dogleg.h"
#include "trust_region_solver_steihuag_cg.h"
#include <zjucad/matrix/io.h>

namespace zjucad
{
namespace LSCO
{

//accepted parameter,default value
//LSCO_solver_in.threshold_outter,1E-4
//LSCO_solver_in.max_iteration_outter,1000
//LSCO_solver_in.threshold_inner,LSCO_solver_in.threshold_outter
//LSCO_solver_in.max_iteration_inner,LSCO_solver_in.max_iteration_outter
//LSCO_solver_in.independent_constraints,false
//LSCO_solver_in.tau,0.995
//LSCO_solver_in.kkt_history_length,5
//LSCO_solver_in.kappa,0.95
//LSCO_solver_in.sigma,0.2
//LSCO_solver_in.rho,0.1
//LSCO_solver_in.eta0,0.1
//LSCO_solver_in.eta1,1/3
//LSCO_solver_in.eta2,2/3
//LSCO_solver_in.tau1,2/3
//LSCO_solver_in.tau2,3/2
//LSCO_solver_in.residual_type,NORM_0
//LSCO_solver_in.initial_mu,1.0
//LSCO_solver_in.adaptive_mu,false	//this is safer but may take more iterations
//LSCO_solver_in.report,false

//returned parameter
//LSCO_solver_out.kkt_violation,value_type
//LSCO_solver_out.kkt_violation_inner,value_type
//LSCO_solver_out.kkt_perturbation,value_type
//LSCO_solver_out.termination_type,std::string
//	successful
//	fail
//LSCO_solver_out.termination_type_detailed,std::string
//	outter_max_iteration
//	inner_max_iteration
//	normal_step_fail
//	predicate_multiplier_fail
//	qp_step_fail
//	outter_threshold_reach
//	minimal_radius_reach

template <typename KERNEL_TYPE>
class LSCO_solver_interior_cg
{
protected:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_kernel_type;
    struct interior_point_hessian : public krylov_matrix<KERNEL_TYPE> {
    public:
        virtual size_t rows() const {
            return m_H_lagrangian->rows()+KERNEL_TYPE::rows(m_ZS);
        }
        virtual void mul(const vector_type& x,vector_type& Ax) {
            m_H_lagrangian->mul(x,Ax);
            size_t nr=m_H_lagrangian->rows();
            size_t nr_inequal=KERNEL_TYPE::rows(m_ZS);
            for(size_t i=0; i<nr_inequal; i++) {
                value_type coef=KERNEL_TYPE::get(i,m_ZS);
                value_type x_val=KERNEL_TYPE::get(nr+i,x);
                KERNEL_TYPE::set(nr+i,coef*x_val,Ax);
            }
        }
        boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_H_lagrangian;
        vector_type m_ZS;
    };
public:
    LSCO_solver_interior_cg(boost::property_tree::ptree& pt)
        :m_opt(pt),
         m_sol_normal(pt.add_child("LSCO_solver_normal",boost::property_tree::ptree())),
         m_sol_qp(pt.add_child("LSCO_solver_qp",boost::property_tree::ptree())) {}
    virtual void solve(objective_function<KERNEL_TYPE>& f, vector_type& x) {
        //variable check
        {
            m_nr=f.nr_inputs();
            m_nr_equal=f.nr_equals();
            m_nr_inequal=f.nr_inequals();
            if(m_nr <= 0 || m_nr_equal < 0 || m_nr_inequal < 0)
                throw "non-positive parameter size";
            if(m_nr != KERNEL_TYPE::rows(x))
                throw "incorrect x size";
        }

        //reset memory
        {
            m_A.reset(new sparse_matrix_type);
            m_rsA.reset(new sparse_matrix_type);
            m_H.reset(new interior_point_hessian);
            KERNEL_TYPE::resize(m_nr_inequal,m_s);
            KERNEL_TYPE::resize(m_nr_inequal,((interior_point_hessian*)m_H.get())->m_ZS);
            KERNEL_TYPE::resize(m_nr+m_nr_inequal,m_rhs);
            KERNEL_TYPE::resize(m_nr+m_nr_inequal,m_new_x);
            KERNEL_TYPE::resize(m_nr+m_nr_inequal,m_g);
            KERNEL_TYPE::resize(m_nr_equal+m_nr_inequal,m_neg_yz);
            KERNEL_TYPE::resize(m_nr_equal+m_nr_inequal,m_A_rhs);
            KERNEL_TYPE::resize(m_nr_equal+m_nr_inequal,m_res_cons);
            KERNEL_TYPE::resize(m_nr_equal+m_nr_inequal,m_nr+m_nr_inequal,*m_A);
        }

        //reset parameters mu
        m_mu=m_opt.get<value_type>("LSCO_solver_in.initial_mu",1.0f);
        bool is_free_mode=m_opt.get<bool>("LSCO_solver_in.adaptive_mu",true);
        std::vector<value_type> kkt_history;

        //reset trust region radius
        init_trust_region_radius(f,x);

        //reset inequality constraint S
        f.eval_constraint_value(x,m_res_cons,true);
        check_objective_consistency(m_res_cons,m_nr_equal+m_nr_inequal);
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type s_val=KERNEL_TYPE::get(m_nr_equal+i,m_res_cons);
            KERNEL_TYPE::set(i,std::max<value_type>(1E-3f,s_val),m_s);
        }

        //reset termination condition
        m_vio_global=std::numeric_limits<value_type>::max();
        m_vio_local=std::numeric_limits<value_type>::max();
        m_vio_local_nrm2=std::numeric_limits<value_type>::max();

        //main outter loop
        value_type thres=m_opt.get<value_type>("LSCO_solver_in.threshold_outter",1E-4f);
        size_t max_iter=m_opt.get<size_t>("LSCO_solver_in.max_iteration_outter",1000);
        //possible report
        if(m_opt.get<bool>("LSCO_solver_in.report",false))
            std::cout << "-------------------------------------------Primal/Interior-Point/CG Begin" << std::endl;
        bool updated=false;
        for(m_iter=0; m_iter < max_iter && m_vio_global > thres; m_iter++) {
            //possible report
            if(m_opt.get<bool>("LSCO_solver_in.report",false))
                std::cout << "outter iter: " << m_iter << " kkt_violation: " << m_vio_global << " mu: " << m_mu << std::endl;
            value_type thres_inner=m_opt.get<value_type>("LSCO_solver_in.threshold_inner",thres);
            size_t max_iter_inner=m_opt.get<size_t>("LSCO_solver_in.max_iteration_inner",max_iter);
            m_first_time=true;
            for(m_iter_inner=0; m_iter_inner < max_iter_inner; m_iter_inner++) {
                //possible report
                if(updated && m_opt.get<bool>("LSCO_solver_in.report",false))
                    std::cout << "\tinner iter: " << m_iter_inner << " kkt_violation: " << m_vio_local << " mu: " << m_mu << std::endl;
                //solve
                if(!solve_inner(f,x,updated,thres_inner))
                    return;
                if(m_vio_local < thres_inner)
                    break;
                //test merit function for trust radius update
                updated=trust_region_update(f,x);
                if(m_radius < m_radius_min) {
                    set_output("successful","minimal_radius_reach");
                    return;
                }
                //update mu adaptively
                if(updated)
                    update_mu_adaptively(is_free_mode,kkt_history);
            }
            if(m_iter_inner == max_iter_inner) {
                set_output("fail","inner_max_iteration");
                return;
            }
            //update mu globally
            update_mu_globally(is_free_mode,kkt_history);
        }
        if(m_iter == max_iter)
            set_output("fail","outter_max_iteration");
        else set_output("successful","outter_threshold_reach");
    }
    void report() const {
        std::cout << "kkt_violation: " << m_opt.get<value_type>("LSCO_solver_out.kkt_violation") << std::endl;
        std::cout << "kkt_violation_inner: " << m_opt.get<value_type>("LSCO_solver_out.kkt_violation_inner") << std::endl;
        std::cout << "kkt_perturbation: " << m_opt.get<value_type>("LSCO_solver_out.kkt_perturbation") << std::endl;
        std::cout << "termination_type: " << m_opt.get<std::string>("LSCO_solver_out.termination_type") << std::endl;
        std::cout << "detailed: " << m_opt.get<std::string>("LSCO_solver_out.termination_type_detailed") << std::endl;
    }
protected:
    //substep
    //input: none, use: m_rhs,m_A_rhs,m_neg_yz, result: m_new_x,m_res_cons,m_A_rhs
    bool solve_normal_step(objective_function<KERNEL_TYPE>& f,vector_type& x) {
        //build A
        KERNEL_TYPE::clear(m_builder);
        f.eval_constraint_gradient(x,m_builder,false);
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type s_val=KERNEL_TYPE::get(i,m_s);
            KERNEL_TYPE::add(m_nr_equal+i,m_nr+i,-s_val,m_builder);
        }
        KERNEL_TYPE::set(m_nr_equal+m_nr_inequal,m_nr+m_nr_inequal,m_builder,*m_A);
        KERNEL_TYPE::row_scale_coeff(*m_A,m_A_rhs,1E-9f);
        *m_rsA=*m_A;
        KERNEL_TYPE::scale_row(m_A_rhs,*m_rsA);
        //build rhs
        f.eval_constraint_value(x,m_res_cons,true);
        check_objective_consistency(m_res_cons,m_nr_equal+m_nr_inequal);
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type s_val=KERNEL_TYPE::get(i,m_s);
            value_type c_val=KERNEL_TYPE::get(m_nr_equal+i,m_res_cons);
            KERNEL_TYPE::set(m_nr_equal+i,c_val-s_val,m_res_cons);
        }
        //solve
        if(m_nr_equal+m_nr_inequal > 0) {
            //build solver
            bool non_singular=m_opt.get<bool>("LSCO_solver_in.independent_constraints",false);
            m_opt.put<bool>("LSCO_solver_normal.trust_region_solver_in.nonsingular_preconditioner",non_singular);
            m_sol_normal.set_hessian_ATA(*m_rsA);
            //apply positive barrier constraint
			value_type tau=m_opt.get<value_type>("LSCO_solver_in.tau",0.995f);
			if(tau > 1.0f || tau <= 0.0f) {
				std::cerr << "invalid tau value" << std::endl;
				return false;
			}
			m_sol_normal.clear_bound();
			for(size_t i=0; i<m_nr_inequal; i++)
				m_sol_normal.set_bound(m_nr+i,-tau*0.5f);		//this parameter is fixed
			try {
                m_sol_normal.assemble();
            } catch(...) {
                return false;
            }
            //solve
            KERNEL_TYPE::copy(m_res_cons,m_neg_yz);
            KERNEL_TYPE::cmul(m_A_rhs,m_neg_yz);
            KERNEL_TYPE::mul_t(*m_rsA,m_neg_yz,m_rhs);
            m_sol_normal.solve(m_rhs,m_radius*0.8f,m_new_x);	//this parameter is fixed
            //check status
            if(m_opt.get<std::string>("LSCO_solver_normal.trust_region_solver_out.termination_type","fail") == "fail")
                return false;
            else return true;
        }
        return true;
    }
    //input: m_new_x, use: m_A_rhs, result: m_neg_yz,m_rhs
    bool predicate_multiplier(objective_function<KERNEL_TYPE>& f,vector_type& x) {
        //build rhs
        f.eval_func_gradient(x,m_rhs,true);
        check_objective_consistency(m_rhs,m_nr+m_nr_inequal);
        for(size_t i=0; i<m_nr_inequal; i++)
            KERNEL_TYPE::set(m_nr+i,-m_mu,m_rhs);
        //solve
        if(m_nr_equal+m_nr_inequal > 0) {
            bool non_singular=m_opt.get<bool>("LSCO_solver_in.independent_constraints",false);
            try {
                m_solA.reset(new kkt_preconditioner_schur_A<KERNEL_TYPE>(m_opt.add_child("LSCO_solver_multiplier",boost::property_tree::ptree()),*m_rsA,non_singular));
            } catch(...) {
                return false;
            }
            KERNEL_TYPE::mul(*m_rsA,m_rhs,m_solA->m_rhs);
            m_solA->m_sol.solve(m_solA->m_rhs,m_neg_yz);
            KERNEL_TYPE::cmul(m_A_rhs,m_neg_yz);
            //project
            for(size_t i=0; i<m_nr_inequal; i++) {
                value_type s_val=KERNEL_TYPE::get(i,m_s);
                assert(s_val > 0.0f && m_mu > 0.0f);
                value_type z_prj=m_mu/s_val;
                KERNEL_TYPE::set(m_nr_equal+i,z_prj,m_neg_yz);
            }
            KERNEL_TYPE::scal(-1.0f,m_neg_yz);
            return true;
        }
        return true;
    }
    //input: m_neg_yz,m_rhs,m_new_x, use: m_new_x,m_neg_yz,m_rhs, result: m_new_x
    bool solve_qp_step(objective_function<KERNEL_TYPE>& f,vector_type& x) {
        //build H
        interior_point_hessian* H_inner=(interior_point_hessian*)m_H.get();
        f.eval_lagrangian_hessian(x,m_neg_yz,1.0f,H_inner->m_H_lagrangian,true);
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type z_val=-KERNEL_TYPE::get(m_nr_equal+i,m_neg_yz);
            value_type s_val=KERNEL_TYPE::get(i,m_s);
            KERNEL_TYPE::set(i,std::max<value_type>(z_val*s_val,1E-3f),H_inner->m_ZS);
        }
        //build solver
        bool non_singular=m_opt.get<bool>("LSCO_solver_in.independent_constraints",false);
        //transfer data to internal cg solver
        m_opt.put<bool>("LSCO_solver_qp.trust_region_solver_in.warm_start",true);
        m_opt.put<bool>("LSCO_solver_qp.trust_region_solver_in.threshold_relative",true);
        m_opt.put<bool>("LSCO_solver_qp.trust_region_solver_in.nonsingular_preconditioner",non_singular);
        m_opt.put<size_t>("LSCO_solver_qp.trust_region_solver_in.max_iteration",20);
        m_opt.put<value_type>("LSCO_solver_qp.trust_region_solver_in.threshold",1E-2f);
		m_sol_qp.set_hessian(m_H);
        //apply positive barrier constraint
        value_type tau=m_opt.get<value_type>("LSCO_solver_in.tau",0.995f);
        if(tau > 1.0f || tau <= 0.0f) {
            std::cerr << "invalid tau value" << std::endl;
            return false;
        }
		m_sol_qp.clear_bound();
		for(size_t i=0; i<m_nr_inequal; i++)
			m_sol_qp.set_bound(m_nr+i,-tau);
        if(m_nr_equal+m_nr_inequal > 0)
            m_sol_qp.set_preconditioner_inv(m_solA,true);
        try {
            m_sol_qp.assemble();
        } catch(...) {
            return false;
        }
		//solve
        KERNEL_TYPE::copy(m_new_x,m_g);
        m_sol_qp.solve(m_rhs,m_radius,m_new_x);
        //check status
        if(m_opt.get<std::string>("LSCO_solver_qp.trust_region_solver_out.termination_type","fail") == "fail")
            return false;
        //project m_new_x again to avoid numeric instability of conjugate gradient
        if(m_nr_equal+m_nr_inequal > 0) {
            KERNEL_TYPE::axpy(-1.0f,m_g,m_new_x);
            m_solA->mul_add(m_new_x,m_g);
            KERNEL_TYPE::copy(m_g,m_new_x);
        }
        return true;
    }
    virtual bool solve_inner(objective_function<KERNEL_TYPE>& f,vector_type& x,bool updated,value_type thres) {
        //solve normal step
        if(!solve_normal_step(f,x)) {
            set_output("fail","normal_step_fail");
            return false;
        }
        //predicate dual variable
        if(!predicate_multiplier(f,x)) {
            set_output("fail","predicate_multiplier_fail");
            return false;
        }
        //update and check termination condition
        kkt_violation(f,x);
        if(m_vio_local < thres)
            return true;
        //update possible quasi-newton matrix
        if(updated || m_first_time) {
            KERNEL_TYPE::copy(m_rhs,m_g,m_nr);
            KERNEL_TYPE::mul_t_add(*m_A,m_neg_yz,m_g);
            if(m_first_time) {
                f.reset(x,m_g);
				
				//evaluate initial merit function
                KERNEL_TYPE::copy(x,m_g,m_nr);
                for(size_t i=0; i<m_nr_inequal; i++)
                    KERNEL_TYPE::set(m_nr+i,KERNEL_TYPE::get(i,m_s),m_g);
                m_last_merit=eval_merit(f,m_g,true);
            } else f.callback_iteration(x,m_g);
        }
        //solve the trust-region subproblem
        if(!solve_qp_step(f,x)) {
            set_output("fail","qp_step_fail");
            return false;
        }
        m_first_time=false;
        return true;
    }
    //parameter tunning: kkt_perturbation
    virtual void update_mu_adaptively(bool& is_free_mode,std::vector<value_type>& kkt_history) {
        if(!is_free_mode || m_nr_inequal == 0)
            return;

        size_t l_max=m_opt.get<size_t>("LSCO_solver_in.kkt_history_length",5);
        value_type kappa=m_opt.get<value_type>("LSCO_solver_in.kappa",0.95f);

        //check if we can continue the free mode
        if(kkt_history.size() != l_max || m_vio_local_nrm2 < (*std::max_element(kkt_history.begin(),kkt_history.end()))*kappa) {
            //apply LOCO law
            value_type min_sz=std::numeric_limits<value_type>::max();
            value_type sz_total=0.0f;
            for(size_t i=0; i<m_nr_inequal; i++) {
                value_type s_val=KERNEL_TYPE::get(i,m_s);
                value_type z_val=-KERNEL_TYPE::get(m_nr_equal+i,m_neg_yz);
                value_type sz_val=s_val*z_val;
                assert(sz_val > 0.0f);
                if(sz_val < min_sz)min_sz=sz_val;
                sz_total+=sz_val;
            }
            value_type sz_avg=sz_total/(value_type)m_nr_inequal;
            value_type yu=min_sz/sz_avg;
            value_type sigma=0.1f*std::pow(std::min<value_type>(0.05f*(1.0f-yu)/yu,2.0f),3);
            m_mu=sigma*sz_avg;
        } else is_free_mode=false;

        //update kkt_history
        if(kkt_history.size() == l_max)
            kkt_history.erase(kkt_history.begin());
        kkt_history.push_back(m_vio_local_nrm2);
    }
    virtual void update_mu_globally(bool& is_free_mode,const std::vector<value_type>& kkt_history) {
        if(is_free_mode || m_nr_inequal == 0)
            return;

        size_t l_max=m_opt.get<size_t>("LSCO_solver_in.kkt_history_length",5);
        value_type kappa=m_opt.get<value_type>("LSCO_solver_in.kappa",0.95f);

        //check if we can return to free mode
        if(kkt_history.size() == l_max && m_vio_local_nrm2 < (*std::max_element(kkt_history.begin(),kkt_history.end()))*kappa)
            is_free_mode=true;
        else {
            //apply Fiacco-McCormick law
            m_mu*=m_opt.get<value_type>("LSCO_solver_in.sigma",0.2f);
        }
    }
    //parameter tunning: trust_region_radius
    virtual void init_trust_region_radius(objective_function<KERNEL_TYPE>& f,vector_type& x) {
        //evaluate initial trust radius
        KERNEL_TYPE::zero(m_g);
        f.eval_func_gradient(x,m_g,false);
        check_objective_consistency(m_g,m_nr+m_nr_inequal);
        m_radius_max=KERNEL_TYPE::nrm2(m_g);
        m_radius_min=m_radius_max*1E-6f;
		m_radius=m_radius_max*0.1f;

        //initialize constraint penalty
        m_nu=0.0f;
    }
    virtual std::pair<value_type,value_type> eval_merit(objective_function<KERNEL_TYPE>& f,const vector_type& x_s,bool same_x) {
        std::pair<value_type,value_type> merit;
        merit.first=f.eval(x_s,same_x);
        f.eval_constraint_value(x_s,m_A_rhs,same_x);
        check_objective_consistency(m_A_rhs,m_nr_equal+m_nr_inequal);
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type s_val=KERNEL_TYPE::get(m_nr+i,x_s);
            merit.first-=m_mu*std::log(s_val);
            value_type ci_val=KERNEL_TYPE::get(m_nr_equal+i,m_A_rhs)-s_val;
            KERNEL_TYPE::set(m_nr_equal+i,ci_val,m_A_rhs);
        }
        merit.second=KERNEL_TYPE::nrm2(m_A_rhs);
        return merit;
    }
    virtual bool trust_region_update(objective_function<KERNEL_TYPE>& f,vector_type& x) {
        //determine predicated constraint-m reduction
        value_type pmred=KERNEL_TYPE::nrm2(m_res_cons);
        KERNEL_TYPE::mul_add(*m_A,m_new_x,m_res_cons);
        pmred-=KERNEL_TYPE::nrm2(m_res_cons);
		assert(pmred > -1E-3f);
        //determine predicated objective reduction
        m_H->mul(m_new_x,m_g);
        value_type pred=-KERNEL_TYPE::dot(m_g,m_new_x)*0.5f-KERNEL_TYPE::dot(m_rhs,m_new_x);
        //update m_nu
        {
            value_type rho=m_opt.get<value_type>("LSCO_solver_in.rho",0.1f);
            value_type denom=pmred*(rho-1.0f);
            if(std::abs(denom) > 1E-3f && m_nu < pred/denom)
                m_nu=1.1f*pred/denom;
        }
        pred+=m_nu*pmred;
        assert(pred >= -1E-3f);

        //determine actual reduction
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type s_val=KERNEL_TYPE::get(i,m_s);
            value_type ds_val=KERNEL_TYPE::get(m_nr+i,m_new_x);
            KERNEL_TYPE::set(m_nr+i,s_val+s_val*ds_val,m_new_x);
        }
        KERNEL_TYPE::axpy(1.0f,x,m_new_x,m_nr);
        std::pair<value_type,value_type> merit=eval_merit(f,m_new_x,false);
        value_type ared=(m_last_merit.first+m_last_merit.second*m_nu)-(merit.first+merit.second*m_nu);
        if(ared < 0.0f)
            std::cout << "increasing iteration" << std::endl;

        //update radius
        value_type eta0=m_opt.get<value_type>("LSCO_solver_in.eta0",1.0f/10.0f);
        value_type eta1=m_opt.get<value_type>("LSCO_solver_in.eta1",1.0f/3.0f);
        value_type eta2=m_opt.get<value_type>("LSCO_solver_in.eta2",2.0f/3.0f);
        value_type tau1=m_opt.get<value_type>("LSCO_solver_in.tau1",2.0f/3.0f);
        value_type tau2=m_opt.get<value_type>("LSCO_solver_in.tau2",3.0f/2.0f);
        if(ared < 0.0f || ared < eta1*pred)
            m_radius*=tau1;
        else if(ared > eta2*pred){
            m_radius*=tau2;
			m_radius=std::min(m_radius,m_radius_max);
		}

        //determine if <x,m_s> is updated
        bool updated=ared > 0.0f && ared > eta0*pred;
        if(updated) {
            //copy to x and m_s
            KERNEL_TYPE::copy(m_new_x,x,m_nr);
            for(size_t i=0; i<m_nr_inequal; i++) {
                value_type s_val=KERNEL_TYPE::get(m_nr+i,m_new_x);
                KERNEL_TYPE::set(i,s_val,m_s);
            }
            m_last_merit=merit;
        }
        return updated;
    }
    //termination check
    void check_objective_consistency(const vector_type& v,size_t nr) const {
        if(KERNEL_TYPE::rows(v) != nr)
            throw "objective function cannot change size of input vector";
    }
    value_type residual(const vector_type& A,const vector_type& B,int type) const {
        switch(type) {
        case NORM_2:
            return std::max(KERNEL_TYPE::nrm2(A),KERNEL_TYPE::nrm2(B));
        case NORM_1:
            return std::max(KERNEL_TYPE::asum(A),KERNEL_TYPE::asum(B));
        case NORM_0:
            return std::max(KERNEL_TYPE::amax(A),KERNEL_TYPE::amax(B));
        default:
            return std::max(KERNEL_TYPE::nrm2(A),KERNEL_TYPE::nrm2(B));
        }
    }
    virtual void kkt_violation(objective_function<KERNEL_TYPE>& f,const vector_type& x) {
        //build the primal residual
        KERNEL_TYPE::copy(m_rhs,m_g);
        KERNEL_TYPE::mul_t_add(*m_A,m_neg_yz,m_g);

        //build the dual residual
        f.eval_constraint_value(x,m_A_rhs,true);
        check_objective_consistency(m_A_rhs,m_nr_equal+m_nr_inequal);
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type c_val=KERNEL_TYPE::get(m_nr_equal+i,m_A_rhs);
            value_type s_val=KERNEL_TYPE::get(i,m_s);
            KERNEL_TYPE::set(m_nr_equal+i,c_val-s_val,m_A_rhs);
        }
        m_vio_local=residual(m_g,m_A_rhs,m_opt.get<int>("LSCO_solver_in.residual_type",NORM_0));
        m_vio_local_nrm2=residual(m_g,m_A_rhs,NORM_2);

        //build the global error
        for(size_t i=0; i<m_nr_inequal; i++) {
            value_type sz_val=KERNEL_TYPE::get(m_nr+i,m_g)-m_mu;
            KERNEL_TYPE::set(m_nr+i,sz_val,m_g);
        }
        m_vio_global=residual(m_g,m_A_rhs,m_opt.get<int>("LSCO_solver_in.residual_type",NORM_0));
    }
    virtual void set_output(const std::string& type,const std::string& type_detailed) {
        m_opt.put<value_type>("LSCO_solver_out.kkt_violation",m_vio_global);
        m_opt.put<value_type>("LSCO_solver_out.kkt_violation_inner",m_vio_local);
        m_opt.put<value_type>("LSCO_solver_out.kkt_perturbation",m_mu);
        m_opt.put<std::string>("LSCO_solver_out.termination_type",type);
        m_opt.put<std::string>("LSCO_solver_out.termination_type_detailed",type_detailed);
        if(m_opt.get<bool>("LSCO_solver_in.report",false))
            report();
    }
protected:
    //parameter
    boost::property_tree::ptree& m_opt;
    size_t m_nr,m_nr_equal,m_nr_inequal,m_iter,m_iter_inner;
    value_type m_radius_max,m_radius_min,m_radius,m_mu,m_nu;
    bool m_first_time;
    //temporary
    std::pair<value_type,value_type> m_last_merit;
    value_type m_vio_local,m_vio_global,m_vio_local_nrm2;
    //data
    vector_type m_s,m_neg_yz,m_res_cons,m_rhs,m_A_rhs,m_new_x,m_g;
    boost::shared_ptr<sparse_matrix_type> m_A,m_rsA;	//constraint matrix JCE 0 \\ JCI -S
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_H;	//hessian of the lagrangian function
    typename KERNEL_TYPE::matrix_builder m_builder;
    //solver
    boost::shared_ptr<kkt_preconditioner_schur_A<KERNEL_TYPE> > m_solA;
    trust_region_solver_double_dogleg<KERNEL_TYPE> m_sol_normal;
    trust_region_solver_steihaug_cg<KERNEL_TYPE> m_sol_qp;
};

}
}

#endif