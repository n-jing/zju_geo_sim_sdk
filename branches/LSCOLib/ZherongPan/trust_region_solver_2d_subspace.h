#ifndef _TRUST_REGION_SOLVER_2D_SUBSPACE_H_
#define _TRUST_REGION_SOLVER_2D_SUBSPACE_H_

#include "trust_region_solver_double_dogleg.h"

extern "C" {
    int solve_quartic(double c[5], double s[4]);
}
namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename MATRIX_TYPE>
struct direct_solver_kernel;

template <typename KERNEL_TYPE>
class trust_region_solver_2d_subspace : public trust_region_solver_double_dogleg<KERNEL_TYPE>
{
    typedef trust_region_solver_double_dogleg<KERNEL_TYPE> BASE_CLASS;
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_kernel_type;
    //functions
    trust_region_solver_2d_subspace(boost::property_tree::ptree& opt)
        :trust_region_solver_double_dogleg<KERNEL_TYPE>(opt) {}
    virtual void assemble() {
		if(!m_bounds.empty())
			throw "2d subspace solver don't support variable bound";
        trust_region_solver_double_dogleg<KERNEL_TYPE>::assemble();
        KERNEL_TYPE::resize(KERNEL_TYPE::rows(*m_H),m_gcu);
        KERNEL_TYPE::resize(KERNEL_TYPE::rows(*m_H),m_gnu);
    }
    virtual void solve(const vector_type& g,value_type delta,vector_type& x) {
        //safe guard
        m_opt.get_child("trust_region_solver_out").clear();
        m_opt.template put<std::string>("trust_region_solver_out.termination_type","successful");
        m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","global_minima");
        if(KERNEL_TYPE::nrm2(g) < 1E-9f) {
            KERNEL_TYPE::zero(x);
            std::cout << "already minimized" << std::endl;
            return;
        }
        m_opt.template put<std::string>("trust_region_solver_out.termination_type_detailed","on_boundary");

        //solve cauchy and newton point
        solve_cauchy(g,delta,m_gc);
        KERNEL_TYPE::copy(g,m_dir);
        KERNEL_TYPE::scal(-1.0f,m_dir);
        m_sol->mul(m_dir,m_gn);

        //unit orthogonalize
        KERNEL_TYPE::copy(m_gc,m_gcu);
        KERNEL_TYPE::scal(1.0f/KERNEL_TYPE::nrm2(m_gc),m_gcu);
        KERNEL_TYPE::copy(m_gn,m_gnu);
        KERNEL_TYPE::axpy(-KERNEL_TYPE::dot(m_gcu,m_gnu),m_gcu,m_gnu);
        value_type nrm_gn=KERNEL_TYPE::nrm2(m_gnu);
        //the two direction are nearly parallel, revert to cauchy point
        if(std::abs(nrm_gn) < 1E-6f) {
            KERNEL_TYPE::copy(m_gc,x);
            return;
        }
        KERNEL_TYPE::scal(1.0/nrm_gn,m_gnu);
        assert(std::abs(KERNEL_TYPE::dot(m_gcu,m_gnu)) < 1E-6f);

        //build reduced linear system
        KERNEL_TYPE::mul(*m_H,m_gcu,m_dir);
        value_type H11=KERNEL_TYPE::dot(m_gcu,m_dir);
        KERNEL_TYPE::mul(*m_H,m_gnu,m_dir);
        value_type H22=KERNEL_TYPE::dot(m_gnu,m_dir);
        value_type H12=KERNEL_TYPE::dot(m_gcu,m_dir);
        value_type STg[2] = {
            KERNEL_TYPE::dot(m_gcu,g),
            KERNEL_TYPE::dot(m_gnu,g)
        };

        //check the positive definiteness
        //if true, just use double dogleg
        //otherwise use quartic equation
        {
            value_type lambda[2];
            value_type coeff_eigen[3]= {
                H11*H22-H12*H12,
                -H11-H22,
                1.0f
            };
            bool indefinite=false;
            int nr=solve_quadric(coeff_eigen,lambda);
            for(size_t i=0; i<nr; i++)
                if(lambda[i] < 0.0f) {
                    indefinite=true;
                    break;
                }
            if(!indefinite) {
                trust_region_solver_double_dogleg<KERNEL_TYPE>::solve(g,delta,x);
                return;
            }
        }

        //solve indefinite
        {
            //start from cauchy point
            KERNEL_TYPE::copy(m_gc,x);
            value_type min_f=eval(g,x);

            //build quartic equation coef
            value_type coeff_quartic[5];
            calc_coef(H11,H12,H22,delta,STg[0],STg[1],coeff_quartic);

            //do search on boundary
            value_type min_s=std::numeric_limits<value_type>::max();
            value_type s[4];
            int nr=solve_quartic(coeff_quartic,s);
            bool found=false;
            for(int i=0; i<nr; ++i) {
                //safe guard
                if(s[i] > -delta+1E-9f && s[i] < delta-1E-9f) {
                    value_type f_b=eval_inner(s[i],STg[0],STg[1],H11,H12,H22,delta);
                    if(f_b < min_f) {
                        min_f=f_b;
                        min_s=s[i];
                        found=true;
                    }
                }
            }
            if(found) {
                // x = u0 * gp_unit_orth + u1 *pn_unit_orth
                KERNEL_TYPE::zero(x);
                KERNEL_TYPE::axpy(min_s,m_gcu,x);
                KERNEL_TYPE::axpy(std::sqrt(delta*delta-min_s*min_s),m_gnu,x);
            }
        }
        if(m_opt.template get<bool>("trust_region_solver_in.report",false))
            report(g,x);
    }
private:
    static void calc_coef(value_type H11,value_type H12,value_type H22,value_type delta,value_type g1,value_type g2,value_type out[5]) {
        /*
        //user variables needed in the method
        value_type H11;
        value_type H12;
        value_type H22;
        value_type delta;
        value_type g1;
        value_type g2;
        */

        //declare temporary variables
        value_type tt1;
        value_type tt2;
        value_type tt3;
        value_type tt4;
        value_type tt5;

        //calculate temporary variables
        tt1=(delta*delta);
        tt2=(g1*g1);
        tt3=(H12*H12);
        tt4=(H11*H11);
        tt5=(H22*H22);

        //outputs
        out[0]=tt1*tt2-(delta*delta*delta*delta)*tt3;
        out[1]=-2*tt1*g1*H22+2*tt1*g2*H12+2*tt1*g1*H11;
        out[2]=tt1*tt5-2*tt1*H11*H22+4*tt1*tt3+tt1*tt4-(g2*g2)-tt2;
        out[3]=2*g1*H22-4*g2*H12-2*g1*H11;
        out[4]=-tt5+2*H11*H22-4*tt3-tt4;
    }
    value_type eval_inner(const value_type alpha,
                          const value_type STg1,const value_type STg2,
                          const value_type H11,const value_type H12,
                          const value_type H22,const value_type delta) {
        value_type beta=std::sqrt(delta*delta-alpha*alpha);
        return STg1*alpha+STg2*beta+	//linear
               0.5f*alpha*alpha*H11+		//quadratic
               0.5f*beta*beta*H22+
               H12*alpha*beta;
    }
    vector_type m_gcu,m_gnu;
};

}
}


#endif
