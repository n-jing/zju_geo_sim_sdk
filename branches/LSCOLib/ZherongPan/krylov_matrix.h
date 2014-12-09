#ifndef _KRYLOV_MATRIX_H_
#define _KRYLOV_MATRIX_H_

#include <boost/shared_ptr.hpp>

namespace zjucad
{
namespace LSCO
{

template <typename VEC_TYPE>
struct kernel_traits;

template <typename KERNEL_TYPE>
struct krylov_matrix {
    typedef typename KERNEL_TYPE::vector_type vector_type;
    virtual size_t rows() const=0;
    virtual size_t cols() const {
        return rows();
    }
    virtual void mul(const vector_type& x,vector_type& Ax)=0;
};

template <typename KERNEL_TYPE>
struct krylov_matrix_AAT : public krylov_matrix<KERNEL_TYPE> {
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    krylov_matrix_AAT(const sparse_matrix_type& A):m_A(A) {}
    virtual size_t rows() const {
        return KERNEL_TYPE::rows(m_A);
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        vector_type tmp;
        KERNEL_TYPE::resize(KERNEL_TYPE::cols(m_A),tmp);
        KERNEL_TYPE::mulT(m_A,x,tmp);
        KERNEL_TYPE::mul(m_A,tmp,Ax);
    }
    const sparse_matrix_type& m_A;
};

template <typename KERNEL_TYPE>
struct krylov_matrix_ATA : public krylov_matrix<KERNEL_TYPE> {
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    krylov_matrix_ATA(const sparse_matrix_type& A):m_A(A) {}
    virtual size_t rows() const {
        return KERNEL_TYPE::cols(m_A);
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        vector_type tmp;
        KERNEL_TYPE::resize(KERNEL_TYPE::rows(m_A),tmp);
        KERNEL_TYPE::mul(m_A,x,tmp);
        KERNEL_TYPE::mul_t(m_A,tmp,Ax);
    }
    const sparse_matrix_type& m_A;
};

template <typename KERNEL_TYPE>
struct default_krylov_matrix : public krylov_matrix<KERNEL_TYPE> {
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    default_krylov_matrix(const sparse_matrix_type& A):m_A(A) {}
    virtual size_t rows() const {
        return KERNEL_TYPE::rows(m_A);
    }
    virtual size_t cols() const {
        return KERNEL_TYPE::cols(m_A);
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        KERNEL_TYPE::mul(m_A,x,Ax);
    }
    const sparse_matrix_type& m_A;
};

template <typename KERNEL_TYPE>
struct default_krylov_matrix_mem : public krylov_matrix<KERNEL_TYPE> {
    typedef typename KERNEL_TYPE::vector_type vector_type;
    typedef typename KERNEL_TYPE::sparse_matrix_type sparse_matrix_type;
    default_krylov_matrix_mem(boost::shared_ptr<sparse_matrix_type> A):m_A(A) {}
    virtual size_t rows() const {
        return KERNEL_TYPE::rows(*m_A);
    }
    virtual size_t cols() const {
        return KERNEL_TYPE::cols(*m_A);
    }
    virtual void mul(const vector_type& x,vector_type& Ax) {
        KERNEL_TYPE::mul(*m_A,x,Ax);
    }
    boost::shared_ptr<sparse_matrix_type> m_A;
};

}
}

#endif
