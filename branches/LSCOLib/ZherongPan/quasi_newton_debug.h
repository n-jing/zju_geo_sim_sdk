#ifndef _QUASI_NEWTON_DEBUG_H_
#define _QUASI_NEWTON_DEBUG_H_

#include "objective_function.h"
#include "kkt_preconditioner.h"

namespace zjucad
{
namespace LSCO
{

template <typename KERNEL_TYPE,typename PARENT>
class random_quadratic : public PARENT
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //functions
    random_quadratic() {
        KERNEL_TYPE::resize(nr_inputs(),m_g);
        KERNEL_TYPE::resize(nr_inputs(),m_Hx);
        KERNEL_TYPE::resize(nr_inputs(),nr_inputs(),m_H);
        KERNEL_TYPE::resize(nr_inputs(),nr_inputs(),m_H0);
        for(size_t r=0; r<nr_inputs(); r++) {
            KERNEL_TYPE::set(r,rand()/(value_type)RAND_MAX,m_g);
            for(size_t c=0; c<nr_inputs(); c++)
                if(r >= c) {
                    value_type v=rand()/(value_type)RAND_MAX;
                    KERNEL_TYPE::set(r,c,v,m_H);
                    KERNEL_TYPE::set(c,r,v,m_H);

                    v=rand()/(value_type)RAND_MAX;
                    KERNEL_TYPE::set(r,c,v,m_H0);
                    KERNEL_TYPE::set(c,r,v,m_H0);
                }
        }
    }
    virtual size_t nr_inputs() const {
        return 50;
    }
    virtual void eval_func_gradient(const vector_type& x,vector_type& gradient,bool same_x) {
        KERNEL_TYPE::mul(m_H,x,m_Hx);
        KERNEL_TYPE::copy(m_Hx,gradient);
        KERNEL_TYPE::axpy(1.0f,m_g,gradient);
    }
    virtual void eval_func_hessian_0(const vector_type& x,boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& hessian,bool same_x) {
        hessian.reset(new default_krylov_matrix<KERNEL_TYPE>(m_H0));
    }
    //data
    vector_type m_g,m_Hx;
    sparse_matrix_type m_H,m_H0;
};

template <typename KERNEL_TYPE>
class krylov_matrix_inv_BFGS : public krylov_matrix<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //functions
    krylov_matrix_inv_BFGS
    (const std::vector<vector_type> xs,
     const std::vector<vector_type>& gs,value_type delta,size_t row)
        :m_delta(delta),m_row(row) {
        vector_type tmp;
        KERNEL_TYPE::resize(rows(),tmp);
        for(size_t i=1; i<xs.size(); i++) {
            KERNEL_TYPE::copy(xs[i],tmp);
            KERNEL_TYPE::axpy(-1.0f,xs[i-1],tmp);
            m_ss.push_back(tmp);

            KERNEL_TYPE::copy(gs[i],tmp);
            KERNEL_TYPE::axpy(-1.0f,gs[i-1],tmp);
            m_ys.push_back(tmp);
        }
    }
    virtual size_t rows() const {
        return m_row;
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        KERNEL_TYPE::copy(x,Ax);
        std::vector<value_type> alphas;
        for(size_t i=0; i<m_ss.size(); i++) {
            size_t j=m_ss.size()-1-i;
            value_type alpha=1.0f/KERNEL_TYPE::dot(m_ss[j],m_ys[j]);
            alpha*=KERNEL_TYPE::dot(m_ss[j],Ax);
            KERNEL_TYPE::axpy(-alpha,m_ys[j],Ax);
            alphas.push_back(alpha);
        }
        KERNEL_TYPE::scal(1.0f/m_delta,Ax);
        for(size_t i=0; i<m_ss.size(); i++) {
            size_t j=m_ss.size()-1-i;
            value_type beta=1.0f/KERNEL_TYPE::dot(m_ss[i],m_ys[i]);
            beta*=KERNEL_TYPE::dot(m_ys[i],Ax);
            KERNEL_TYPE::axpy(alphas[j]-beta,m_ss[i],Ax);
        }
    }
    //data
    std::vector<vector_type> m_ss,m_ys;
    value_type m_delta;
    size_t m_row;
};

template <typename KERNEL_TYPE>
class krylov_matrix_inv_SR1 : public krylov_matrix_inv_BFGS<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;

    typedef krylov_matrix_inv_BFGS<KERNEL_TYPE> BASE_CLASS;
    //functions
    krylov_matrix_inv_SR1
    (const std::vector<vector_type> xs,
     const std::vector<vector_type>& gs,
     const boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& H0,size_t row)
        :krylov_matrix_inv_BFGS<KERNEL_TYPE>(xs,gs,1.0f,row),m_H0(H0) {
        vector_type Hy;
        KERNEL_TYPE::resize(rows(),Hy);
        for(size_t i=0; i<BASE_CLASS::m_ss.size(); i++) {
            mul(BASE_CLASS::m_ys[i],Hy);
            KERNEL_TYPE::scal(-1.0f,Hy);
            KERNEL_TYPE::axpy(1.0f,BASE_CLASS::m_ss[i],Hy);
            m_sSHy.push_back(Hy);
            m_coef.push_back(1.0f/KERNEL_TYPE::dot(Hy,BASE_CLASS::m_ys[i]));
        }
    }
    virtual size_t rows() const {
        return BASE_CLASS::m_row;
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        m_H0->mul(x,Ax);
        for(size_t i=0; i<m_coef.size(); i++) {
            value_type coef=KERNEL_TYPE::dot(m_sSHy[i],x)*m_coef[i];
            KERNEL_TYPE::axpy(coef,m_sSHy[i],Ax);
        }
    }
    //data
    const boost::shared_ptr<krylov_matrix<KERNEL_TYPE> >& m_H0;
    std::vector<vector_type> m_sSHy;
    std::vector<value_type> m_coef;
};

template <typename KERNEL_TYPE>
void BFGS_debug(objective_function_L_D_BFGS<KERNEL_TYPE>& f,size_t nrc,size_t nr,bool damp)
{
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;

    //these settings are for BFGS
    f.set_nr_correct(nrc);
    f.set_damp(damp);
    std::vector<vector_type> xs,gs;
    for(size_t i=0; i<nr; i++) {
        vector_type new_x,new_g;
        KERNEL_TYPE::resize(f.nr_inputs(),new_x);
        KERNEL_TYPE::resize(f.nr_inputs(),new_g);
        for(size_t id=0; id<f.nr_inputs(); id++)
            KERNEL_TYPE::set(id,rand()/(value_type)RAND_MAX,new_x);

        f.eval_func_gradient(new_x,new_g,false);
        xs.push_back(new_x);
        gs.push_back(new_g);
    }

    //build BFGS matrix
    f.reset(xs[0],gs[0]);
    for(size_t i=1; i<nr; i++) {
        f.callback_iteration(xs[i],gs[i]);
    }

    //limit memory for inverse hessian
    while(xs.size() > f.nr_correct()+1) {
        xs.erase(xs.begin());
        gs.erase(gs.begin());
    }

    //test it with inverse BFGS
    krylov_matrix_inv_BFGS<KERNEL_TYPE> H(xs,gs,f.delta(),f.nr_inputs());

    //test 100 samples
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > hess;
    {
        vector_type verbose;
        f.eval_lagrangian_hessian(verbose,verbose,1.0f,hess,false);
    }
    for(size_t i=0; i<100; i++) {
        vector_type a,Ha,BHa;
        KERNEL_TYPE::resize(f.nr_inputs(),a);
        KERNEL_TYPE::resize(f.nr_inputs(),Ha);
        KERNEL_TYPE::resize(f.nr_inputs(),BHa);
        for(size_t id=0; id<f.nr_inputs(); id++)
            KERNEL_TYPE::set(id,rand()/(value_type)RAND_MAX,a);
        H.mul(a,Ha);
        hess->mul(Ha,BHa);
        KERNEL_TYPE::axpy(-1.0f,a,BHa);
        std::cout << "Err: " << KERNEL_TYPE::nrm2(BHa) << std::endl;
    }
}

template <typename KERNEL_TYPE>
void SR1_debug(random_quadratic<KERNEL_TYPE,objective_function_L_SR1<KERNEL_TYPE> >& f,size_t nrc,size_t nr)
{
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;

    //these settings are for BFGS
    f.set_nr_correct(nrc);
    std::vector<vector_type> xs,gs;
    for(size_t i=0; i<nr; i++) {
        vector_type new_x,new_g;
        KERNEL_TYPE::resize(f.nr_inputs(),new_x);
        KERNEL_TYPE::resize(f.nr_inputs(),new_g);
        for(size_t id=0; id<f.nr_inputs(); id++)
            KERNEL_TYPE::set(id,rand()/(value_type)RAND_MAX,new_x);

        f.eval_func_gradient(new_x,new_g,false);
        xs.push_back(new_x);
        gs.push_back(new_g);
    }

    //build SR1 matrix
    f.reset(xs[0],gs[0]);
    for(size_t i=1; i<nr; i++) {
        f.callback_iteration(xs[i],gs[i]);
    }

    //limit memory for inverse hessian
    while(xs.size() > f.nr_correct()+1) {
        xs.erase(xs.begin());
        gs.erase(gs.begin());
    }

    //test it with inverse SR1
    boost::property_tree::ptree pt;
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > H0(new kkt_preconditioner_M<KERNEL_TYPE>(pt,f.m_H0,false));
    krylov_matrix_inv_SR1<KERNEL_TYPE> H(xs,gs,H0,f.nr_inputs());

    //test 100 samples
    boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > hess;
    {
        vector_type verbose;
        f.eval_lagrangian_hessian(verbose,verbose,1.0f,hess,false);
    }
    for(size_t i=0; i<100; i++) {
        vector_type a,Ha,BHa;
        KERNEL_TYPE::resize(f.nr_inputs(),a);
        KERNEL_TYPE::resize(f.nr_inputs(),Ha);
        KERNEL_TYPE::resize(f.nr_inputs(),BHa);
        for(size_t id=0; id<f.nr_inputs(); id++)
            KERNEL_TYPE::set(id,rand()/(value_type)RAND_MAX,a);
        H.mul(a,Ha);
        hess->mul(Ha,BHa);
        KERNEL_TYPE::axpy(-1.0f,a,BHa);
        std::cout << "Err: " << KERNEL_TYPE::nrm2(BHa) << std::endl;
    }
}

template <typename KERNEL_TYPE>
void quasi_newton_debug(bool BFGS,bool SR1)
{
    if(BFGS) {
        srand(1000);
        random_quadratic<KERNEL_TYPE,objective_function_L_D_BFGS<KERNEL_TYPE> > f;
        srand(1000);
        std::cout << "-------------------------------------------BFGS implementation" << std::endl;
        BFGS_debug(f,40,30,false);

        srand(1000);
        std::cout << "-------------------------------------------L-BFGS implementation" << std::endl;
        BFGS_debug(f,15,30,false);

        srand(1000);
        std::cout << "-------------------------------------------L-D-BFGS implementation" << std::endl;
        BFGS_debug(f,15,30,true);
    }
    if(SR1) {
        srand(1000);
        random_quadratic<KERNEL_TYPE,objective_function_L_SR1<KERNEL_TYPE> > f;
        srand(1000);
        std::cout << "-------------------------------------------SR1 implementation" << std::endl;
        SR1_debug(f,40,30);

        srand(1000);
        std::cout << "-------------------------------------------L-SR1 implementation" << std::endl;
        SR1_debug(f,15,30);
    }
}

}
}

#endif
