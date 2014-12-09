#ifndef _KERNEL_DEBUG_H_
#define _KERNEL_DEBUG_H_

#include "kernel.h"
#include "kernel_zjucad_matrix.h"
#include "kernel_common_file_pzr.h"

#include "krylov_matrix.h"
#include <iostream>
#include <assert.h>

namespace zjucad
{
namespace LSCO
{

template <typename KERNEL_TYPE>
static void kernel_vector_debug()
{
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;

    vector_type av,bv,axpbv,consv,zerov;
    size_t N=50,actual=40;

    KERNEL_TYPE::resize(N,av);
    KERNEL_TYPE::resize(N,bv);
    KERNEL_TYPE::resize(N,axpbv);
    for(size_t i=0; i<actual; i++) {
        KERNEL_TYPE::set(i,rand()*2.0f/(value_type)RAND_MAX-1.0f,av);
        KERNEL_TYPE::set(i,rand()*2.0f/(value_type)RAND_MAX-1.0f,bv);
    }
    for(size_t i=actual; i<N; i++) {
        KERNEL_TYPE::set(i,0.0f,av);
        KERNEL_TYPE::set(i,1.0f,bv);
    }
    std::cout << "dot: " << KERNEL_TYPE::dot(av,bv,actual) << std::endl;
    KERNEL_TYPE::cmul(bv,av,actual);
    std::cout << "cmul-nrm2: " << KERNEL_TYPE::nrm2(av,actual) << std::endl;
    KERNEL_TYPE::scal(2.0f,av,actual);
    std::cout << "dot-after-scale: " << KERNEL_TYPE::dot(av,bv,actual) << std::endl;
    assert(KERNEL_TYPE::rows(av) == N);

    KERNEL_TYPE::copy(bv,axpbv,actual);
    KERNEL_TYPE::axpy(rand()/(value_type)RAND_MAX,av,axpbv,actual);
    std::cout << "axpy-nrm2: " << KERNEL_TYPE::nrm2(axpbv,actual) << std::endl;
    std::cout << "axpy-asum: " << KERNEL_TYPE::asum(axpbv,actual) << std::endl;
    std::cout << "axpy-amax: " << KERNEL_TYPE::amax(axpbv,actual) << std::endl;

    std::cout << "before: " << KERNEL_TYPE::nrm2(av,actual) << " " << KERNEL_TYPE::nrm2(bv,actual) << std::endl;
    KERNEL_TYPE::swap(av,bv,actual);
    std::cout << "after: " << KERNEL_TYPE::nrm2(av,actual) << " " << KERNEL_TYPE::nrm2(bv,actual) << std::endl;

    KERNEL_TYPE::resize(N,consv);
    KERNEL_TYPE::resize(N,zerov);
    KERNEL_TYPE::zero(zerov);
    KERNEL_TYPE::set(rand()*2.0f/(value_type)RAND_MAX-1.0f,consv,actual);
    std::cout << "before-add-zero: " << KERNEL_TYPE::nrm2(av,actual) << std::endl;
    KERNEL_TYPE::axpy(1.0f,zerov,av,actual);
    std::cout << "after-add-zero: " << KERNEL_TYPE::nrm2(av,actual) << std::endl;
    KERNEL_TYPE::axpy(1.0f,consv,av,actual);
    std::cout << "after-add-constant: " << KERNEL_TYPE::nrm2(av,actual) << std::endl;
}

template <typename KERNEL_TYPE>
static void kernel_matrix_debug(bool use_builder)
{
    typedef typename KERNEL_TYPE::value_type value_type;
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type matrix_type;

    vector_type rv,cv,Acv;
    size_t R=50,C=25;
    KERNEL_TYPE::resize(R,rv);
    KERNEL_TYPE::resize(C,cv);
    for(size_t i=0; i<R; i++)
        KERNEL_TYPE::set(i,rand()*2.0f/(value_type)RAND_MAX-1.0f,rv);
    for(size_t i=0; i<C; i++)
        KERNEL_TYPE::set(i,rand()*2.0f/(value_type)RAND_MAX-1.0f,cv);
    std::cout << "arbitrary-pos-query: " << KERNEL_TYPE::get(rand()%KERNEL_TYPE::rows(rv),rv) << std::endl;

    matrix_type A;
    if(use_builder) {
        KERNEL_TYPE::matrix_builder builder;
        KERNEL_TYPE::clear(builder);
        for(size_t r=0; r<R; r++)
            for(size_t c=0; c<C; c++) {
                size_t nn=rand()%3;
                for(size_t n=0; n<nn; n++)
                    KERNEL_TYPE::add(r,c,rand()*2.0f/(value_type)RAND_MAX-1.0f,builder);
            }
        KERNEL_TYPE::set(R,C,builder,A);
        KERNEL_TYPE::set(R,C,builder,A);
    } else {
        KERNEL_TYPE::resize(R,C,A);
        for(size_t r=0; r<R; r++)
            for(size_t c=0; c<C; c++)
                KERNEL_TYPE::set(r,c,rand()*2.0f/(value_type)RAND_MAX-1.0f,A);
    }
    assert(KERNEL_TYPE::rows(A) == R);
    assert(KERNEL_TYPE::cols(A) == C);

    KERNEL_TYPE::resize(R,Acv);
    KERNEL_TYPE::mul(A,cv,Acv);
    std::cout << "Acv-norm1: " << KERNEL_TYPE::asum(Acv) << std::endl;

    vector_type row_scale;
    KERNEL_TYPE::resize(R,row_scale);
    KERNEL_TYPE::row_scale_coeff(A,row_scale,1E-9f);
    KERNEL_TYPE::scale_row(row_scale,A);
    KERNEL_TYPE::mul(A,cv,Acv);
    std::cout << "Acv-norm1: " << KERNEL_TYPE::asum(Acv) << std::endl;

    default_krylov_matrix<KERNEL_TYPE> dkm(A);
    assert(dkm.rows() == R);
    dkm.mul(cv,Acv);
    std::cout << "Acv-norm1-krylov: " << KERNEL_TYPE::asum(Acv) << std::endl;

    KERNEL_TYPE::resize(R,Acv);
    KERNEL_TYPE::set(rand()*2.0f/(value_type)RAND_MAX-1.0f,Acv);
    KERNEL_TYPE::mul_add(A,cv,Acv);
    std::cout << "Acv-add-norm1: " << KERNEL_TYPE::asum(Acv) << std::endl;

    KERNEL_TYPE::resize(R,Acv);
    KERNEL_TYPE::set(rand()*2.0f/(value_type)RAND_MAX-1.0f,Acv);
    KERNEL_TYPE::mul_sub(A,cv,Acv);
    std::cout << "Acv-sub-norm1: " << KERNEL_TYPE::asum(Acv) << std::endl;

    KERNEL_TYPE::resize(C,Acv);
    KERNEL_TYPE::mul_t(A,rv,Acv);
    std::cout << "ATcv-norm1: " << KERNEL_TYPE::asum(Acv) << std::endl;

    KERNEL_TYPE::resize(C,Acv);
    KERNEL_TYPE::set(rand()*2.0f/(value_type)RAND_MAX-1.0f,Acv);
    KERNEL_TYPE::mul_t_add(A,rv,Acv);
    std::cout << "ATcv-add-norm1: " << KERNEL_TYPE::asum(Acv) << std::endl;

    KERNEL_TYPE::resize(C,Acv);
    KERNEL_TYPE::set(rand()*2.0f/(value_type)RAND_MAX-1.0f,Acv);
    KERNEL_TYPE::mul_t_sub(A,rv,Acv);
    std::cout << "ATcv-sub-norm1: " << KERNEL_TYPE::asum(Acv) << std::endl;
}

template <typename SCALAR_TYPE>
void kernel_debug()
{
    srand(1000);
    typedef kernel_traits<plain_matrix<SCALAR_TYPE> > KT1;
    std::cout << "-------------------------------------------plain implementation" << std::endl;
    kernel_vector_debug<KT1>();
    std::cout << std::endl;
    kernel_matrix_debug<KT1>(true);
    std::cout << std::endl;
    kernel_matrix_debug<KT1>(false);
    std::cout << std::endl;

    srand(1000);
    typedef kernel_traits<hj::sparse::csc<SCALAR_TYPE,typename hj::sparse::idx_type> > KT2;
    std::cout << "-------------------------------------------zjucad::matrix implementation" << std::endl;
    kernel_vector_debug<KT2>();
    std::cout << std::endl;
    kernel_matrix_debug<KT2>(true);
    std::cout << std::endl;
    kernel_matrix_debug<KT2>(false);
    std::cout << std::endl;

    srand(1000);
    typedef kernel_traits<COMMON::FixedSparseMatrix<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > > KT3;
    std::cout << "-------------------------------------------CommonFilePZR implementation" << std::endl;
    kernel_vector_debug<KT3>();
    std::cout << std::endl;
    kernel_matrix_debug<KT3>(true);
    std::cout << std::endl;
    kernel_matrix_debug<KT3>(false);
    std::cout << std::endl;
}

}
}

#endif
