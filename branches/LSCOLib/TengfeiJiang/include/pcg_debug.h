#ifndef _PCG_DEBUG_H_
#define _PCG_DEBUG_H_

#include "kernel_common_file_pzr.h"
#include "direct_solver_kernel_common_file_pzr.h"

#include "pcg_solver.h"
#include "kkt_preconditioner.h"

#include <CommonFilePZR/solvers/LinearSolver.h>	//as reference
#include <CommonFilePZR/solvers/PMinres.h>	//as reference

namespace zjucad
{
namespace LSCO
{

template <typename SCALAR_TYPE>
void pcg_debug(size_t n,size_t m,size_t nrC,size_t prob0)
{
    //the test is specific to common file for it's convenient
    //other implementations should perform same by kernel_debug.h
    typedef kernel_traits<COMMON::FixedSparseMatrix<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > > KERNEL_TYPE;
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type matrix_type;

    matrix_type A,C,invADiag,ADiag;
    vector_type b,result,resultRef;
    A.resize(n*m,n*m);
    C.resize(nrC,n*m);
    invADiag.resize(n*m,n*m);
    ADiag.resize(n*m,n*m);
    b.resize(n*m);
    result.resize(n*m);
    resultRef.resize(n*m);

    //build model problem
    for(sizeType r=0; r<(sizeType)n; r++)
        for(sizeType c=0; c<(sizeType)m; c++) {
#define GI(R,C) ((R)*m+(C))
            if(prob0 == 0) {
                A.addToElement(GI(r,c),GI(r,c),rand()/(value_type)RAND_MAX+1.0f);
                invADiag.addToElement(GI(r,c),GI(r,c),sqrt(1.0f/A(GI(r,c),GI(r,c))));
                ADiag.addToElement(GI(r,c),GI(r,c),sqrt(A(GI(r,c),GI(r,c))));
            } else {
                if(r == 0) {
                    //NEUMANN
                } else {
                    A.addToElement(GI(r,c),GI(r,c),1.0f);
                    A.addToElement(GI(r,c),GI(r-1,c),-1.0f);
                }
                if(r == n-1) {
                    //NEUMANN
                } else {
                    A.addToElement(GI(r,c),GI(r,c),1.0f);
                    A.addToElement(GI(r,c),GI(r+1,c),-1.0f);
                }
                if(c == 0) {
                    //NEUMANN
                } else {
                    A.addToElement(GI(r,c),GI(r,c),1.0f);
                    A.addToElement(GI(r,c),GI(r,c-1),-1.0f);
                }
                if(c == m-1) {
                    //DIRICHLET
                    A.addToElement(GI(r,c),GI(r,c),1.0f);
                } else {
                    A.addToElement(GI(r,c),GI(r,c),1.0f);
                    A.addToElement(GI(r,c),GI(r,c+1),-1.0f);
                }
                invADiag.addToElement(GI(r,c),GI(r,c),1.0f/A(GI(r,c),GI(r,c)));
                ADiag.addToElement(GI(r,c),GI(r,c),sqrt(A(GI(r,c),GI(r,c))));
            }
            b[GI(r,c)]=rand()*2.0f/(value_type)RAND_MAX-1.0f;
#undef GI
        }

    //build constraints
    for(size_t c=0; c<nrC; c++) {
        std::map<size_t,SCALAR_TYPE> cmap;
        for(size_t cv=0; cv<4; cv++)
            cmap[rand()%(n*m)]=rand()*2.0f/(value_type)RAND_MAX-1.0f;
        for(typename std::map<size_t,SCALAR_TYPE>::const_iterator
                beg=cmap.begin(),end=cmap.end(); beg!=end; beg++)
            C.setElement(c,beg->first,beg->second);
    }

    boost::property_tree::ptree pt;
    pt.put<bool>("pcg_in.warm_start",false);
    pt.put<int>("pcg_in.residual_type",NORM_0);
    pt.put<SCALAR_TYPE>("pcg_in.threshold",1E-8f);
    pt.put<bool>("pcg_in.flexible",true);
    pt.put<bool>("pcg_in.report",true);

    //TEST A: based CG test
    std::cout << "-------------------------------------------CG Test" << std::endl;
    //solve using CommonFilePZR
    {
        COMMON::PCGSolver<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE>,COMMON::NoPreconSolver<SCALAR_TYPE> > sol;
        sol.setMatrix(A,true);
        boost::shared_ptr<COMMON::Callback<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > > cb;
        cb.reset(new COMMON::Callback<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> >);
        sol.setCallback(cb);
        sol.solve(b,resultRef);
    }
    //solve using LSCO::pcg_solver
    {
        pcg_solver<KERNEL_TYPE> sol(pt);
        sol.set_A(A);
        sol.solve(b,result);
        sol.report();
    }
    //TEST B: preconditioner test
    std::cout << "-------------------------------------------PCG Test" << std::endl;
    {
        COMMON::PCGSolver<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE>,COMMON::ExplicitPreconSolver<SCALAR_TYPE> > sol;
        sol.setMatrix(A,false);
        sol.getPre()->setMatrix(invADiag,false);
        boost::shared_ptr<COMMON::Callback<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > > cb;
        cb.reset(new COMMON::Callback<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> >);
        sol.setCallback(cb);
        sol.solve(b,resultRef);
    }
    //solve using LSCO::pcg_solver
    {
        pcg_solver<KERNEL_TYPE> sol(pt);
        sol.set_A(A);
        sol.set_pre(invADiag);
        sol.solve(b,result);
        sol.report();
    }
    //TEST C: let's add some constraints
    std::cout << "-------------------------------------------KKT-PCG Test" << std::endl;
    {
        vector_type bAug,resultAug;
        bAug.resize(n*m+nrC);
        bAug.setZero();
        bAug.block(0,0,n*m,1)=b;

        resultAug.resize(n*m+nrC);
        resultAug.setZero();

        COMMON::PMINRESSolver<SCALAR_TYPE> sol;
        boost::shared_ptr<COMMON::KKTKrylovMatrix<SCALAR_TYPE> > krylov(new COMMON::KKTKrylovMatrix<SCALAR_TYPE>(A,C));
        sol.setKrylovMatrix(krylov);
        sol.setSolverParameters(1E-8f,1000);
        sol.solve(bAug,resultAug);
        resultRef=resultAug.block(0,0,n*m,1);
    }
    {
        pcg_solver<KERNEL_TYPE> sol(pt);
        sol.set_A(A);
        boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > precon(new kkt_preconditioner_schur_A<KERNEL_TYPE>(pt,C,false));
        sol.set_pre(precon);
        sol.solve(b,result);
        sol.report();
        std::cout << "Err: " << (resultRef-result).norm() << " of " << result.norm() << std::endl;
    }
    {
        pcg_solver<KERNEL_TYPE> sol(pt);
        sol.set_A(A);
        boost::shared_ptr<krylov_matrix<KERNEL_TYPE> > precon(new kkt_preconditioner_M_A<KERNEL_TYPE>(pt,ADiag,C,false));
        sol.set_pre(precon);
        sol.solve(b,result);
        sol.report();
        std::cout << "Err: " << (resultRef-result).norm() << " of " << result.norm() << std::endl;
    }
}

}
}

#endif
