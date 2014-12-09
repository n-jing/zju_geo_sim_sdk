#ifndef _DIRECT_SOLVER_KERNEL_DEBUG_H_
#define _DIRECT_SOLVER_KERNEL_DEBUG_H_

#include "direct_solver_kernel.h"
#include "direct_solver_kernel_zjucad_matrix.h"
#include "direct_solver_kernel_common_file_pzr.h"

namespace zjucad
{
namespace LSCO
{

template <typename SOL_TYPE>
void solve_debug(const std::string& info,SOL_TYPE& sol)
{
    typedef typename SOL_TYPE::kernel_type kernel_type;
    typedef typename SOL_TYPE::value_type value_type;
    typedef typename SOL_TYPE::vector_type vector_type;
    typedef typename SOL_TYPE::sparse_matrix_type matrix_type;

    vector_type rhs,result;
    kernel_type::resize(sol.rows(),rhs);
    kernel_type::resize(sol.rows(),result);
    for(size_t i=0; i<sol.rows(); i++)
        kernel_type::set(i,rand()*2.0f/(value_type)RAND_MAX-1.0f,rhs);
    sol.solve(rhs,result);
    std::cout << info << " norm2: " << kernel_type::nrm2(result) << std::endl;
}

template <typename DIRECT_SOLVER_TYPE>
void direct_solver_solve_debug()
{
    typedef typename DIRECT_SOLVER_TYPE::kernel_type kernel_type;
    typedef typename DIRECT_SOLVER_TYPE::value_type value_type;
    typedef typename DIRECT_SOLVER_TYPE::vector_type vector_type;
    typedef typename DIRECT_SOLVER_TYPE::sparse_matrix_type matrix_type;

    matrix_type A;
    matrix_type H;
    matrix_type CM;
    size_t N=50,M=25,C=5;
    value_type v;
    boost::property_tree::ptree pt;
    DIRECT_SOLVER_TYPE sol(pt);

    kernel_type::resize(N,M,A);
    for(size_t r=0; r<N; r++)
        for(size_t c=0; c<M; c++)
            kernel_type::set(r,c,rand()*2.0f/(value_type)RAND_MAX-1.0f,A);
    sol.build_AAT(A,false,false);
    direct_solver_kernel_report(pt);
    solve_debug("ATA-nonsingular",sol);

    sol.build_AAT(A,false,true);
    direct_solver_kernel_report(pt);
    solve_debug("ATA-nonsingular-known",sol);

    sol.build_AAT(A,true,false);
    direct_solver_kernel_report(pt);
    solve_debug("AAT-singular",sol);

    kernel_type::resize(N,N,H);
    for(size_t r=0; r<N; r++)
        for(size_t c=0; c<N; c++) {
            if(r<=c) {
                v=rand()*2.0f/(value_type)RAND_MAX-1.0f;
                kernel_type::set(r,c,v,H);
                kernel_type::set(c,r,v,H);
            }
        }
    sol.build_H(H,false,false);
    direct_solver_kernel_report(pt);
    solve_debug("H-nonsingular",sol);

    sol.build_H(H,false,true);
    direct_solver_kernel_report(pt);
    solve_debug("H-nonsingular-known",sol);

    kernel_type::resize(C,N,CM);
    for(size_t r=0; r<C; r++)
        for(size_t c=0; c<N; c++)
            kernel_type::set(r,c,rand()*2.0f/(value_type)RAND_MAX-1.0f,CM);
    sol.build_KKT(H,CM,false);
    direct_solver_kernel_report(pt);
    solve_debug("KKT-singular",sol);

    sol.build_KKT_schur(H,CM,false);
    direct_solver_kernel_report(pt);
    solve_debug("KKT-schur-singular",sol);
}

template <typename SCALAR_TYPE>
static void direct_solver_debug()
{
    srand(1000);
    typedef direct_solver_kernel<plain_matrix<SCALAR_TYPE> > DKT1;
    std::cout << "-------------------------------------------plain implementation" << std::endl;
    direct_solver_solve_debug<DKT1>();

    srand(1000);
    typedef direct_solver_kernel<hj::sparse::csc<SCALAR_TYPE> > DKT2;
    std::cout << "-------------------------------------------zjucad::matrix implementation" << std::endl;
    direct_solver_solve_debug<DKT2>();

    srand(1000);
    typedef direct_solver_kernel<COMMON::FixedSparseMatrix<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > > DKT3;
    std::cout << "-------------------------------------------CommonFilePZR implementation" << std::endl;
    direct_solver_solve_debug<DKT3>();
}

}
}

#endif
