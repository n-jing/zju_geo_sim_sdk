#ifndef TRUST_REGION_SOLVER_DEBUG_H
#define TRUST_REGION_SOLVER_DEBUG_H

#include "kernel_common_file_pzr.h"
#include "direct_solver_kernel_common_file_pzr.h"
#include "trust_region_solver_steihuag_cg.h"
#include "trust_region_solver_double_dogleg.h"
#include "trust_region_solver_2d_subspace.h"
#include <set>

namespace zjucad
{
namespace LSCO
{

template <typename KERNEL_TYPE>
class func
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //interface
    virtual size_t n() const=0;
    virtual value_type eval(const vector_type& x,vector_type& g,sparse_matrix_type& h)=0;
    virtual void eval_A(boost::shared_ptr<sparse_matrix_type>& A) {
        A.reset((sparse_matrix_type*)NULL);
    }
    virtual void init(vector_type& x)=0;
};
template <typename KERNEL_TYPE>
class func_Bohachecsky : public func<KERNEL_TYPE>
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //interface
    virtual size_t n() const {
        return 2;
    }
    virtual value_type eval(const vector_type& x,vector_type& g,sparse_matrix_type& h) {
        value_type x0=KERNEL_TYPE::get(0,x);
        value_type x1=KERNEL_TYPE::get(1,x);
        KERNEL_TYPE::resize(2,g);
        KERNEL_TYPE::resize(2,2,h);
        KERNEL_TYPE::set(0,2.0f*x0+0.9f*M_PI*sin(3.0f*M_PI*x0),g);
        KERNEL_TYPE::set(1,4.0f*x1+1.6f*M_PI*sin(4.0f*M_PI*x1),g);
        KERNEL_TYPE::set(0,0,2.0f+2.7f*M_PI*M_PI*cos(3.0f*M_PI*x0),h);
        KERNEL_TYPE::set(1,1,4.0f+6.4f*M_PI*M_PI*cos(4.0f*M_PI*x1),h);
        return x0*x0+2.0f*x1*x1-0.3f*cos(3.0f*M_PI*x0)-0.4f*cos(4.0f*M_PI*x1)+0.7f;
    }
    virtual void init(vector_type& x) {
        KERNEL_TYPE::resize(2,x);
        KERNEL_TYPE::set(0,-5.0f,x);
        KERNEL_TYPE::set(1,-1.0f,x);
    }
};
template <typename KERNEL_TYPE>
class func_Simplest0 : public func<KERNEL_TYPE>
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //interface
    virtual size_t n() const {
        return 1;
    }
    virtual value_type eval(const vector_type& x,vector_type& g,sparse_matrix_type& h) {
        value_type x0=KERNEL_TYPE::get(0,x);
        KERNEL_TYPE::resize(1,g);
        KERNEL_TYPE::resize(1,1,h);
        KERNEL_TYPE::set(0,4.0f*x0*x0*x0,g);
        KERNEL_TYPE::set(0,0,12.0f*x0*x0,h);
        return x0*x0*x0*x0;
    }
    virtual void init(vector_type& x) {
        KERNEL_TYPE::resize(1,x);
        KERNEL_TYPE::set(0,-0.5f,x);
    }
};
template <typename KERNEL_TYPE>
class func_Simplest1 : public func<KERNEL_TYPE>
{
public:
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    //interface
    virtual size_t n() const {
        return 20;
    }
    virtual value_type eval(const vector_type& x,vector_type& g,sparse_matrix_type& h) {
        value_type ret=0.0f;
        KERNEL_TYPE::resize(20,g);
        KERNEL_TYPE::resize(20,20,h);
        for(size_t i=0; i<20; i++) {
            value_type xi=KERNEL_TYPE::get(i,x);
            value_type entry=(value_type)i+1.0f;
            KERNEL_TYPE::set(i,i,entry,h);
            ret+=xi*xi*entry*0.5f;
            KERNEL_TYPE::set(i,xi*entry,g);
        }
        return ret;
    }
    virtual void init(vector_type& x) {
        KERNEL_TYPE::resize(20,x);
        KERNEL_TYPE::set(1.0f,x);
    }
};

template <typename KERNEL_TYPE>
void trust_region_solve(func<KERNEL_TYPE>& f,trust_region_solver<KERNEL_TYPE>& sol,typename KERNEL_TYPE::value_type delta,size_t n,bool bound)
{
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    vector_type g,x,dx;
    boost::shared_ptr<sparse_matrix_type> A,h(new sparse_matrix_type);

    KERNEL_TYPE::resize(f.n(),g);
    KERNEL_TYPE::resize(f.n(),x);
    KERNEL_TYPE::resize(f.n(),dx);
    f.init(x);
    f.eval_A(A);

    for(size_t i=0; i<n; i++) {
        std::cout << "iter: " << i << std::endl;
        f.eval(x,g,*h);
		if(bound)
		{
			sol.clear_bound();
			std::set<size_t> id;
			for(size_t i=0;i<KERNEL_TYPE::rows(x)/10;i++){
				size_t bd_id=rand()%KERNEL_TYPE::rows(x);
				while(id.find(bd_id) != id.end())
					bd_id=rand()%KERNEL_TYPE::rows(x);
				sol.set_bound(bd_id,-delta/100.0f);
			}
		}
        sol.set_hessian(h);
        sol.set_constraint(A);
        sol.assemble();
        sol.solve(g,delta,dx);
        assert(KERNEL_TYPE::nrm2(dx) < delta+1E-5f);
        KERNEL_TYPE::axpy(1.0f,dx,x);
        std::cout << "func-nonlinear: " << f.eval(x,g,*h) << std::endl;
    }
}

template <typename KERNEL_TYPE>
void trust_region_solver_debug(bool cg,bool dd,bool s2d,bool bound,bool report=false)
{
    boost::property_tree::ptree pt;
    std::cout << "-------------------------------------------Steihaug-CG implementation" << std::endl;
    if(cg) {
        trust_region_solver_steihaug_cg<KERNEL_TYPE> sol_cg(pt);
        pt.put<bool>("trust_region_solver_in.report",report);

        std::cout << "func_Bohachecsky" << std::endl;
        func_Bohachecsky<KERNEL_TYPE> f1;
        trust_region_solve<KERNEL_TYPE>(f1,sol_cg,5.0f,10,bound);
        std::cout << std::endl;

        std::cout << "func_Simplest0" << std::endl;
        func_Simplest0<KERNEL_TYPE> f2;
        trust_region_solve<KERNEL_TYPE>(f2,sol_cg,0.05f,10,bound);
        std::cout << std::endl;

        std::cout << "func_Simplest1" << std::endl;
        func_Simplest1<KERNEL_TYPE> f3;
        trust_region_solve<KERNEL_TYPE>(f3,sol_cg,4.0f,2,bound);
        std::cout << std::endl;

        std::cout << "func_Simplest1" << std::endl;
        trust_region_solve<KERNEL_TYPE>(f3,sol_cg,4.1f,2,bound);
        std::cout << std::endl;
    }
    std::cout << "-------------------------------------------Double Dogleg implementation" << std::endl;
    if(dd) {
        trust_region_solver_double_dogleg<KERNEL_TYPE> sol_dd(pt);
        pt.put<bool>("trust_region_solver_in.report",report);

        std::cout << "func_Simplest0" << std::endl;
        func_Simplest0<KERNEL_TYPE> f2;
        trust_region_solve<KERNEL_TYPE>(f2,sol_dd,0.06f,10,bound);
        std::cout << std::endl;

        std::cout << "func_Simplest1" << std::endl;
        func_Simplest1<KERNEL_TYPE> f3;
        trust_region_solve<KERNEL_TYPE>(f3,sol_dd,3.9f,2,bound);
        std::cout << std::endl;
    }
    std::cout << "-------------------------------------------2d-subspace implementation" << std::endl;
    if(s2d) {
        trust_region_solver_2d_subspace<KERNEL_TYPE> sol_2d(pt);
        pt.put<bool>("trust_region_solver_in.report",report);

        std::cout << "func_Bohachecsky" << std::endl;
        func_Bohachecsky<KERNEL_TYPE> f1;
        trust_region_solve<KERNEL_TYPE>(f1,sol_2d,15.0f,10,bound);
        std::cout << std::endl;

        std::cout << "func_Simplest0" << std::endl;
        func_Simplest0<KERNEL_TYPE> f2;
        trust_region_solve<KERNEL_TYPE>(f2,sol_2d,0.05f,10,bound);
        std::cout << std::endl;

        std::cout << "func_Simplest1" << std::endl;
        func_Simplest1<KERNEL_TYPE> f3;
        trust_region_solve<KERNEL_TYPE>(f3,sol_2d,3.9f,2,bound);
        std::cout << std::endl;
    }
}

}
}

#endif
