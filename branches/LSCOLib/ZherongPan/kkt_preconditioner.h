#ifndef _KKT_PRECONDITIOENR_H_
#define _KKT_PRECONDITIOENR_H_

#include "krylov_matrix.h"
#include "direct_solver_kernel_report.h"
#include <boost/property_tree/ptree.hpp>

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename MATRIX_TYPE>
struct direct_solver_kernel;

template <typename KERNEL_TYPE>
class kkt_preconditioner_schur_A : public krylov_matrix<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_type;
    //functions
    kkt_preconditioner_schur_A(boost::property_tree::ptree& pt,const sparse_matrix_type& A,bool non_singular)
        :m_sol(pt.add_child("kkt_preconditioner_schur_A",boost::property_tree::ptree())),m_A(A) {
        m_sol.build_AAT(A,true,non_singular);
        if(pt.get<std::string>("kkt_preconditioner_schur_A.direct_solver_out.termination_type","fail") == "fail")
            throw "direct solver build fail";
        //direct_solver_kernel_report(pt.get_child("kkt_preconditioner_schur_A"));
        KERNEL_TYPE::resize(KERNEL_TYPE::rows(A),m_rhs);
        KERNEL_TYPE::resize(KERNEL_TYPE::rows(A),m_lambda);
    }
    virtual size_t rows() const {
        return KERNEL_TYPE::cols(m_A);
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        KERNEL_TYPE::mul(m_A,x,m_rhs);
        m_sol.solve(m_rhs,m_lambda);
        KERNEL_TYPE::copy(x,Ax);
        KERNEL_TYPE::mul_t_sub(m_A,m_lambda,Ax);
    }
    virtual void mul_add(const vector_type& x,vector_type& Ax) {
        KERNEL_TYPE::mul(m_A,x,m_rhs);
        m_sol.solve(m_rhs,m_lambda);
        KERNEL_TYPE::axpy(1.0f,x,Ax);
        KERNEL_TYPE::mul_t_sub(m_A,m_lambda,Ax);
    }
    //data
    vector_type m_rhs,m_lambda;
    const sparse_matrix_type& m_A;
    direct_solver_type m_sol;
};

template <typename KERNEL_TYPE>
class kkt_preconditioner_ATA : public krylov_matrix<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_type;
    //functions
    kkt_preconditioner_ATA(boost::property_tree::ptree& pt,const sparse_matrix_type& A,bool non_singular)
        :m_sol(pt.add_child("kkt_preconditioner_ATA",boost::property_tree::ptree())) {
        m_sol.build_AAT(A,false,non_singular);
        if(pt.get<std::string>("kkt_preconditioner_ATA.direct_solver_out.termination_type","fail") == "fail")
            throw "direct solver build fail";
        //direct_solver_kernel_report(pt.get_child("kkt_preconditioner_ATA"));
    }
    virtual size_t rows() const {
        return m_sol.rows();
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        m_sol.solve(x,Ax);
    }
    //data
    direct_solver_type m_sol;
};

template <typename KERNEL_TYPE>
class kkt_preconditioner_M_A : public krylov_matrix<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_type;
    //functions
    kkt_preconditioner_M_A(boost::property_tree::ptree& pt,const sparse_matrix_type& M,const sparse_matrix_type& A,bool non_singular)
        :m_sol(pt.add_child("kkt_preconditioner_M_A",boost::property_tree::ptree())),m_M(M),m_A(A) {
        m_sol.build_KKT(M,A,non_singular);
        if(pt.get<std::string>("kkt_preconditioner_M_A.direct_solver_out.termination_type","fail") == "fail")
            throw "direct solver build fail";
        //direct_solver_kernel_report(pt.get_child("kkt_preconditioner_M_A"));
        KERNEL_TYPE::resize(rows()+KERNEL_TYPE::rows(m_A),m_rhs);
        KERNEL_TYPE::resize(rows()+KERNEL_TYPE::rows(m_A),m_x);
        KERNEL_TYPE::zero(m_rhs);
    }
    virtual size_t rows() const {
        return KERNEL_TYPE::rows(m_M);
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        KERNEL_TYPE::copy(x,m_rhs,rows());
        m_sol.solve(m_rhs,m_x);
        KERNEL_TYPE::copy(m_x,Ax,rows());
    }
    //data
    vector_type m_rhs,m_x;
    const sparse_matrix_type& m_M;
    const sparse_matrix_type& m_A;
    direct_solver_type m_sol;
};

template <typename KERNEL_TYPE>
class kkt_preconditioner_invM_A : public krylov_matrix<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_type;
    //functions
    kkt_preconditioner_invM_A(boost::property_tree::ptree& pt,const sparse_matrix_type& invM,const sparse_matrix_type& A,bool non_singular)
        :m_sol(pt.add_child("kkt_preconditioner_invM_A",boost::property_tree::ptree())),m_invM(invM),m_A(A) {
        m_sol.build_KKT_schur(invM,A,non_singular);
        if(pt.get<std::string>("kkt_preconditioner_invM_A.direct_solver_out.termination_type","fail") == "fail")
            throw "direct solver build fail";
        //direct_solver_kernel_report(pt.get_child("kkt_preconditioner_invM_A"));
        KERNEL_TYPE::resize(KERNEL_TYPE::rows(A),m_rhs);
        KERNEL_TYPE::resize(KERNEL_TYPE::rows(A),m_lambda);
        KERNEL_TYPE::resize(KERNEL_TYPE::cols(A),m_tmp);
    }
    virtual size_t rows() const {
        return KERNEL_TYPE::rows(m_invM);
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        KERNEL_TYPE::mul(m_invM,x,Ax);
        KERNEL_TYPE::mul(m_A,Ax,m_rhs);
        m_sol.solve(m_rhs,m_lambda);

        KERNEL_TYPE::copy(x,m_tmp);
        KERNEL_TYPE::mul_t_sub(m_A,m_lambda,m_tmp);
        KERNEL_TYPE::mul(m_invM,m_tmp,Ax);
    }
    //data
    vector_type m_rhs,m_lambda,m_tmp;
    const sparse_matrix_type& m_invM;
    const sparse_matrix_type& m_A;
    direct_solver_type m_sol;
};

template <typename KERNEL_TYPE>
class kkt_preconditioner_M : public krylov_matrix<KERNEL_TYPE>
{
public:
    //typedefs
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    typedef direct_solver_kernel<sparse_matrix_type> direct_solver_type;
    //functions
    kkt_preconditioner_M(boost::property_tree::ptree& pt,const sparse_matrix_type& M,bool non_singular)
        :m_sol(pt.add_child("kkt_preconditioner_M",boost::property_tree::ptree())) {
        m_sol.build_H(M,non_singular,true);
        if(pt.get<std::string>("kkt_preconditioner_M.direct_solver_out.termination_type","fail") == "fail")
            throw "direct solver build fail";
        //direct_solver_kernel_report(pt.get_child("kkt_preconditioner_M"));
    }
    virtual size_t rows() const {
        return m_sol.rows();
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        m_sol.solve(x,Ax);
    }
    //data
    direct_solver_type m_sol;
};

}
}

#endif
