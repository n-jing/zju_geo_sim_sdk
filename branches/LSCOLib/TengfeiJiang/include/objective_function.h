#ifndef _OBJECTIVE_FUNCTION_H_
#define _OBJECTIVE_FUNCTION_H_

#include "krylov_matrix.h"
#include "direct_solver_dense_kernel.h"

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename MATRIX_TYPE>
struct direct_solver_kernel;

template <typename KERNEL_TYPE>
class objective_function
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //objective functions
    objective_function():m_hessian(new default_krylov_matrix_mem<KERNEL_TYPE>(boost::shared_ptr<sparse_matrix_type>())) {}
    virtual size_t nr_inputs() const {
        assert(false);
        return 0;
    }
    virtual size_t nr_equals() const {
        return 0;
    }
    virtual size_t nr_inequals() const {
        return 0;
    }
    virtual value_type eval(const vector_type& x,bool same_x) {
        assert(false);
        return 0.0f;
    }
    virtual void eval_func_gradient(const vector_type& x,vector_type& gradient,bool same_x) {
        assert(false);
    }
    virtual void eval_lagrangian_hessian(const vector_type& x,const vector_type& lambda,value_type coef,boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& hessian,bool same_x) {
        hessian=m_hessian;
        if(!m_hessian->m_A)
            m_hessian->m_A.reset(new sparse_matrix_type);
        KERNEL_TYPE::clear(m_builder);
        eval_lagrangian_hessian(x,lambda,coef,m_builder,same_x);
        KERNEL_TYPE::set(nr_inputs(),nr_inputs(),m_builder,*(m_hessian->m_A));
    }
    virtual void eval_lagrangian_hessian(const vector_type& x,const vector_type& lambda,value_type coef,typename KERNEL_TYPE::matrix_builder& hessian,bool same_x) {
        assert(false);
    }
    virtual void eval_constraint_value(const vector_type& x,vector_type& c,bool same_x) {}
    virtual void eval_constraint_gradient(const vector_type& x,typename KERNEL_TYPE::matrix_builder& c_gradient,bool same_x) {}
    //internal function, don't call these
    virtual void reset(const vector_type& x,const vector_type& gradient) {}
    virtual void callback_iteration(const vector_type& x,const vector_type& gradient) {}
protected:
    boost::shared_ptr<default_krylov_matrix_mem<KERNEL_TYPE> > m_hessian;
    typename KERNEL_TYPE::matrix_builder m_builder;
};

template <typename KERNEL_TYPE>
class objective_function_L_D_BFGS : public objective_function<KERNEL_TYPE>
{
protected:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_dense_kernel<value_type> dense_solver_type;

    typedef objective_function<KERNEL_TYPE> BASE_CLASS;
    //updation formula (internal function, don't call these)
    class quasi_newton_krylov_matrix_L_D_BFGS : public krylov_matrix<KERNEL_TYPE>
    {
    public:
        quasi_newton_krylov_matrix_L_D_BFGS
        (std::vector<boost::shared_ptr<vector_type> >& S,
         std::vector<boost::shared_ptr<vector_type> >& Y,size_t rows)
            :m_S(S),m_Y(Y),m_rows(rows) {}
        virtual size_t rows() const {
            return m_rows;
        }
        virtual void mul(const vector_type& x,vector_type& Ax) {
            KERNEL_TYPE::copy(x,Ax);
            KERNEL_TYPE::scal(m_delta,Ax);

            size_t nr=m_S.size();
            for(size_t i=0; i<nr; i++) {
                m_sol.set_rhs(i,m_delta*KERNEL_TYPE::dot(*(m_S[i]),x));
                m_sol.set_rhs(i+nr,KERNEL_TYPE::dot(*(m_Y[i]),x));
            }
            m_sol.solve();
            for(size_t i=0; i<nr; i++) {
                KERNEL_TYPE::axpy(-m_delta*m_sol.get_result(i),*(m_S[i]),Ax);
                KERNEL_TYPE::axpy(-m_sol.get_result(i+nr),*(m_Y[i]),Ax);
            }
        }
        virtual bool build() {
            value_type v;
            size_t nr=m_S.size();
            m_sol.reset(nr*2);
            for(size_t r=0; r<nr; r++)
                for(size_t c=0; c<nr; c++) {
                    if(r >= c) {
                        v=m_delta*KERNEL_TYPE::dot(*(m_S[r]),*(m_S[c]));
                        if(r == c) {
                            //delta*S^TS
                            m_sol.set(r,r,v);
                            //-D
                            v=-KERNEL_TYPE::dot(*(m_S[r]),*(m_Y[r]));
                            m_sol.set(r+nr,r+nr,v);
                        } else {
                            //delta*S^TS
                            m_sol.set(r,c,v);
                            m_sol.set(c,r,v);
                            //L
                            v=KERNEL_TYPE::dot(*(m_S[r]),*(m_Y[c]));
                            m_sol.set(r,c+nr,v);
                            m_sol.set(c+nr,r,v);
                        }
                    }
                }
            return m_sol.build(true,false);
        }
        //data
        std::vector<boost::shared_ptr<vector_type> >& m_S;
        std::vector<boost::shared_ptr<vector_type> >& m_Y;
        dense_solver_type m_sol;
        value_type m_delta;
        size_t m_rows;
    };
public:
    objective_function_L_D_BFGS() {
        set_damp(true);
    }
    size_t nr_correct() const {
        return m_correct;
    }
    value_type delta() const {
        return m_hessian_qn->m_delta;
    }
    void set_nr_correct(size_t correct) {
        assert(correct > 0 && correct <=nr_inputs());
        m_correct=correct;
    }
    void set_damp(bool damp) {
        m_damp=damp;
    }
    virtual void eval_lagrangian_hessian(const vector_type& x,const vector_type& lambda,value_type coef,boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& hessian,bool same_x) {
        hessian=m_hessian_qn;
    }
    virtual void reset(const vector_type& x,const vector_type& gradient) {
        m_hessian_qn.reset(new quasi_newton_krylov_matrix_L_D_BFGS(m_S,m_Y,nr_inputs()));
        m_hessian_qn->m_delta=1.0f;
        KERNEL_TYPE::resize(nr_inputs(),m_lx);
        KERNEL_TYPE::resize(nr_inputs(),m_lg);
        KERNEL_TYPE::copy(x,m_lx);
        KERNEL_TYPE::copy(gradient,m_lg);
        KERNEL_TYPE::resize(nr_inputs(),m_s);
        KERNEL_TYPE::resize(nr_inputs(),m_y);
        m_S.clear();
        m_Y.clear();
    }
    virtual void callback_iteration(const vector_type& x,const vector_type& gradient) {
        callback_iteration_inner(x,gradient);
        KERNEL_TYPE::copy(x,m_lx);
        KERNEL_TYPE::copy(gradient,m_lg);
    }
protected:
    virtual void callback_iteration_inner(const vector_type& x,const vector_type& gradient) {
        //calculate m_s and m_y
        KERNEL_TYPE::copy(x,m_s,nr_inputs());
        KERNEL_TYPE::axpy(-1.0f,m_lx,m_s);
        KERNEL_TYPE::copy(gradient,m_y,nr_inputs());
        KERNEL_TYPE::axpy(-1.0f,m_lg,m_y);

        //safe-guard
        value_type ys=KERNEL_TYPE::dot(m_s,m_y);
        if(std::abs(ys) < 1E-6f) {
            //update last x and gradient
            std::cout << "skipping update" << std::endl;
            return;
        }

        //limit memory
        if(m_S.size() == m_correct) {
            m_S.erase(m_S.begin());
            m_Y.erase(m_Y.begin());
        }

        //damp m_y if required (see: "numerical optimization" proc(18.2))
        if(m_damp) {
            vector_type& bs=m_lx;	//reused
            if(!m_hessian_qn->build())
                throw "error building quasi-newton matrix";
            m_hessian_qn->mul(m_s,bs);
            value_type sbs=KERNEL_TYPE::dot(bs,m_s);
            assert(sbs > 0.0f);	//this is the merit we obtained by damping
            if(ys < 0.2f*sbs) {
                value_type theta=0.8f*sbs/(sbs-ys);
                KERNEL_TYPE::scal(theta,m_y);
                KERNEL_TYPE::axpy(1.0f-theta,bs,m_y);
            }
        }
        m_hessian_qn->m_delta=
            KERNEL_TYPE::dot(m_y,m_y)/
            KERNEL_TYPE::dot(m_s,m_y);

        //build compact formula
        m_S.push_back(boost::shared_ptr<vector_type>(new vector_type));
        m_Y.push_back(boost::shared_ptr<vector_type>(new vector_type));
        KERNEL_TYPE::resize(nr_inputs(),*(m_S.back()));
        KERNEL_TYPE::resize(nr_inputs(),*(m_Y.back()));
        KERNEL_TYPE::copy(m_s,*(m_S.back()));
        KERNEL_TYPE::copy(m_y,*(m_Y.back()));
        if(!m_hessian_qn->build())
            throw "error building quasi-newton matrix";
    }
    //data: compact representation L-BFGS matrix
    bool m_damp;
    size_t m_correct;
    vector_type m_lx,m_lg,m_s,m_y;
    std::vector<boost::shared_ptr<vector_type> > m_S;
    std::vector<boost::shared_ptr<vector_type> > m_Y;
    boost::shared_ptr<quasi_newton_krylov_matrix_L_D_BFGS> m_hessian_qn;
};

template <typename KERNEL_TYPE>
class objective_function_L_SR1 : public objective_function_L_D_BFGS<KERNEL_TYPE>
{
protected:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_dense_kernel<value_type> dense_solver_type;

    typedef objective_function_L_D_BFGS<KERNEL_TYPE> BASE_CLASS;
    //updation formula
    class quasi_newton_krylov_matrix_L_SR1 : public BASE_CLASS::quasi_newton_krylov_matrix_L_D_BFGS
    {
    public:
        typedef typename BASE_CLASS::quasi_newton_krylov_matrix_L_D_BFGS BASE_QUASI_NEWTON;
        quasi_newton_krylov_matrix_L_SR1
        (std::vector<boost::shared_ptr<vector_type> >& S,
         std::vector<boost::shared_ptr<vector_type> >& Y,
         boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& h0,size_t rows)
            :BASE_QUASI_NEWTON(S,Y,rows),m_hessian_0(h0) {}
        virtual size_t rows() const {
            return BASE_QUASI_NEWTON::m_rows;
        }
        virtual void mul(const vector_type& x,vector_type& Ax) {
            if(m_hessian_0)
                m_hessian_0->mul(x,Ax);
            else KERNEL_TYPE::copy(x,Ax);
            size_t nr=BASE_QUASI_NEWTON::m_S.size();
            for(size_t i=0; i<nr; i++)
                BASE_QUASI_NEWTON::m_sol.set_rhs(i,KERNEL_TYPE::dot(*(BASE_QUASI_NEWTON::m_Y[i]),x));
            BASE_QUASI_NEWTON::m_sol.solve();
            for(size_t i=0; i<nr; i++)
                KERNEL_TYPE::axpy(BASE_QUASI_NEWTON::m_sol.get_result(i),*(BASE_QUASI_NEWTON::m_Y[i]),Ax);
        }
        virtual bool build() {
            value_type v;
            size_t nr=BASE_QUASI_NEWTON::m_S.size();
            BASE_QUASI_NEWTON::m_sol.reset(nr);
            for(size_t r=0; r<nr; r++)
                for(size_t c=0; c<nr; c++) {
                    if(r > c) {
                        v=KERNEL_TYPE::dot(*(BASE_QUASI_NEWTON::m_S[r]),*(BASE_QUASI_NEWTON::m_Y[c]));
                        BASE_QUASI_NEWTON::m_sol.set(r,c,v);
                        BASE_QUASI_NEWTON::m_sol.set(c,r,v);
                    } else if(r == c) {
                        v=KERNEL_TYPE::dot(*(BASE_QUASI_NEWTON::m_S[r]),*(BASE_QUASI_NEWTON::m_Y[r]));
                        BASE_QUASI_NEWTON::m_sol.set(r,r,v);
                    }
                }
            return BASE_QUASI_NEWTON::m_sol.build(true,false);
        }
        //data
        boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& m_hessian_0;
    };
public:
    virtual void reset(const vector_type& x,const vector_type& gradient) {
        objective_function_L_D_BFGS<KERNEL_TYPE>::reset(x,gradient);
        BASE_CLASS::m_hessian_qn.reset(new quasi_newton_krylov_matrix_L_SR1(BASE_CLASS::m_S,BASE_CLASS::m_Y,m_hessian_0,nr_inputs()));
        eval_func_hessian_0(x,m_hessian_0,false);
    }
    virtual void eval_func_hessian_0(const vector_type& x,boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& hessian,bool same_x) {
        hessian.reset((krylov_matrix<KERNEL_TYPE>*)NULL);
    }
protected:
    virtual void callback_iteration_inner(const vector_type& x,const vector_type& gradient) {
        //calculate m_s and m_y
        KERNEL_TYPE::copy(x,BASE_CLASS::m_s,nr_inputs());
        KERNEL_TYPE::axpy(-1.0f,BASE_CLASS::m_lx,BASE_CLASS::m_s);
        KERNEL_TYPE::copy(gradient,BASE_CLASS::m_y,nr_inputs());
        KERNEL_TYPE::axpy(-1.0f,BASE_CLASS::m_lg,BASE_CLASS::m_y);

        //limit memory
        boost::shared_ptr<vector_type> m_S0;
        boost::shared_ptr<vector_type> m_Y0;
        if(BASE_CLASS::m_S.size() == BASE_CLASS::m_correct) {
            m_S0=BASE_CLASS::m_S[0];
            m_Y0=BASE_CLASS::m_Y[0];
            BASE_CLASS::m_S.erase(BASE_CLASS::m_S.begin());
            BASE_CLASS::m_Y.erase(BASE_CLASS::m_Y.begin());
        }

        //safe-guard to avoid unboundedness (see: "numerical optimization" eqn(6.26))
        vector_type& bs=BASE_CLASS::m_lx;	//reused
        value_type R=1E-6f;
        {
            if(m_S0) {
                if(!BASE_CLASS::m_hessian_qn->build())
                    throw "error building quasi-newton matrix";
            }
            BASE_CLASS::m_hessian_qn->mul(BASE_CLASS::m_s,bs);
            KERNEL_TYPE::scal(-1.0f,bs);
            KERNEL_TYPE::axpy(1.0f,BASE_CLASS::m_y,bs);
            if(std::abs(KERNEL_TYPE::dot(BASE_CLASS::m_s,bs)) < R*KERNEL_TYPE::nrm2(BASE_CLASS::m_s)*KERNEL_TYPE::nrm2(bs)) {
                if(m_S0) {
                    BASE_CLASS::m_S.insert(BASE_CLASS::m_S.begin(),m_S0);
                    BASE_CLASS::m_Y.insert(BASE_CLASS::m_Y.begin(),m_Y0);
                    if(!BASE_CLASS::m_hessian_qn->build())
                        throw "error building quasi-newton matrix";
                }
                std::cout << "skipping update" << std::endl;
                return;
            }
        }

        //build compact formula
        BASE_CLASS::m_S.push_back(boost::shared_ptr<vector_type>(new vector_type));
        BASE_CLASS::m_Y.push_back(boost::shared_ptr<vector_type>(new vector_type));
        KERNEL_TYPE::resize(nr_inputs(),*(BASE_CLASS::m_S.back()));
        KERNEL_TYPE::resize(nr_inputs(),*(BASE_CLASS::m_Y.back()));
        KERNEL_TYPE::copy(BASE_CLASS::m_s,*(BASE_CLASS::m_S.back()));
        KERNEL_TYPE::copy(BASE_CLASS::m_y,*(BASE_CLASS::m_Y.back()));
        if(m_hessian_0)
            m_hessian_0->mul(BASE_CLASS::m_s,bs);
        else KERNEL_TYPE::copy(BASE_CLASS::m_s,bs);
        KERNEL_TYPE::axpy(-1.0f,bs,*(BASE_CLASS::m_Y.back()));

        eval_func_hessian_0(x,m_hessian_0,false);
        if(!BASE_CLASS::m_hessian_qn->build())
            throw "error building quasi-newton matrix";
    }
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > m_hessian_0;
};

}
}

#endif
