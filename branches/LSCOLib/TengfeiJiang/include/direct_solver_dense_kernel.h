#ifndef _DIRECT_SOLVER_DENSE_KERNEL_H_
#define _DIRECT_SOLVER_DENSE_KERNEL_H_

#include "kernel.h"
#include "direct_solver_kernel.h"

namespace zjucad
{
namespace LSCO
{

template <typename SCALAR_TYPE>
struct direct_solver_dense_kernel {
public:
    //build system
    typedef kernel_traits<plain_matrix<SCALAR_TYPE> > kernel_type;
    direct_solver_dense_kernel():m_sol(m_opt) {}
    void reset(size_t n) {
        kernel_type::resize(n,n,m_A);
        kernel_type::resize(n,m_rhs);
        kernel_type::resize(n,m_result);
    }
    void set(size_t r,size_t c,SCALAR_TYPE a) {
        kernel_type::set(r,c,a,m_A);
    }
    void set_rhs(size_t r,SCALAR_TYPE a) {
        kernel_type::set(r,a,m_rhs);
    }
    //solve
    bool build(bool non_singular,bool positive) {
        m_sol.build_H(m_A,non_singular,positive);
        return m_opt.get<std::string>("direct_solver_out.termination_type","fail") == "successful";
    }
    void solve() {
        m_sol.solve(m_rhs,m_result);
    }
    SCALAR_TYPE get_result(size_t r) {
        return kernel_type::get(r,m_result);
    }
protected:
    boost::property_tree::ptree m_opt;
    direct_solver_kernel<plain_matrix<SCALAR_TYPE> > m_sol;
    plain_matrix<SCALAR_TYPE> m_A;
    plain_vector<SCALAR_TYPE> m_rhs,m_result;
};

}
}

#endif