#ifndef _KERNEL_COMMON_FILE_PZR_H_
#define _KERNEL_COMMON_FILE_PZR_H_

#include <stdint.h>
#include <CommonFilePZR/solvers/MatVec.h>

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename SCALAR_TYPE>
struct kernel_traits<COMMON::FixedSparseMatrix<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > > {
    //typedefs
    typedef typename COMMON::Kernel<SCALAR_TYPE> kernel_inner;
    typedef SCALAR_TYPE value_type;
    typedef typename COMMON::Kernel<value_type>::Vec vector_type;
    typedef typename COMMON::FixedSparseMatrix<value_type,kernel_inner> sparse_matrix_type;
    typedef typename std::vector<Eigen::Triplet<value_type,sizeType> > matrix_builder;
    //functions of blas: plain slow implementation
    static value_type dot(const vector_type& x,const vector_type& y,int64_t n=-1) {
        return kernel_inner::dot(x,y,n);
    }
    static void axpy(value_type a,const vector_type& x,vector_type& y,int64_t n=-1) {
        kernel_inner::addScaled(a,x,y,n);
    }
    static void copy(const vector_type& x,vector_type& y,int64_t n=-1) {
        kernel_inner::copy(x,y,n);
    }
    static void swap(vector_type& x,vector_type& y,int64_t n=-1) {
        vector_type tmp=x;
        kernel_inner::copy(y,x,n);
        kernel_inner::copy(tmp,y,n);
    }
    static value_type nrm2(const vector_type& x,int64_t n=-1) {
        return kernel_inner::norm(x,n);
    }
    static value_type asum(const vector_type& x,int64_t n=-1) {
        if(n == -1) {
            n=x.size();
        } else {
            assert((size_t)x.size() > (size_t)n);
        }
        value_type r = 0;
        for(sizeType i=0; i<(sizeType)n; ++i)
            r+=std::abs(x[i]);
        return r;
    }
    static void scal(value_type a,vector_type& x,int64_t n=-1) {
        kernel_inner::scale(a,x,n);
    }
    static value_type amax(const vector_type& x,int64_t n=-1) {
        return kernel_inner::absMax(x,n);
    }
    //other necessary functions
    static void resize(size_t n,vector_type& x) {
        x.resize(n);
    }
    static void zero(vector_type& x,int64_t n=-1) {
        set(0.0f,x,n);
    }
    static void set(value_type a,vector_type& x,int64_t n=-1) {
        kernel_inner::set(x,a,n);
    }
    static void set(size_t i,value_type a,vector_type& x) {
        x[i]=a;
    }
    static value_type get(size_t i,const vector_type& x) {
        return x[i];
    }
    static size_t rows(const vector_type& x) {
        return x.size();
    }
    static void cmul(const vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1)
            y.array()*=x.array();
        else y.block(0,0,n,1).array()*=x.block(0,0,n,1).array();
    }
    //matrix vector multiplication function
    static void mul(const sparse_matrix_type& A,const vector_type& x,vector_type& Ax) {
        A.multiply(x,Ax);
    }
    static void mul_add(const sparse_matrix_type& A,const vector_type& x,vector_type& bPAx) {
        A.multiplyAdd(x,bPAx);
    }
    static void mul_sub(const sparse_matrix_type& A,const vector_type& x,vector_type& bSAx) {
        A.multiplySubtract(x,bSAx);
    }
    static void mul_t(const sparse_matrix_type& A,const vector_type& x,vector_type& ATx) {
        A.multiplyTranspose(x,ATx);
    }
    static void mul_t_add(const sparse_matrix_type& A,const vector_type& x,vector_type& bPATx) {
        A.multiplyTransposeAdd(x,bPATx);
    }
    static void mul_t_sub(const sparse_matrix_type& A,const vector_type& x,vector_type& bSATx) {
        A.multiplyTransposeSubtract(x,bSATx);
    }
    //other necessary matrix functions
    static void resize(size_t r,size_t c,sparse_matrix_type& x) {
        x.resize(r,c);
    }
    static void set(size_t r,size_t c,value_type a,sparse_matrix_type& x) {
        x.setElement(r,c,a);
    }
    static size_t rows(const sparse_matrix_type& x) {
        return x.rows();
    }
    static size_t cols(const sparse_matrix_type& x) {
        return x.cols();
    }
    static void row_scale_coeff(const sparse_matrix_type& x,vector_type& s,value_type thres) {
        const std::vector<sizeType>& off=x.getRowOffset();
        const typename sparse_matrix_type::ROW& row=x.getValue();
        for(sizeType r=0; r<x.rows(); r++) {
#define SQR(x) (x)*(x)
            value_type& coeff=s[r];
            coeff=0.0f;
            for(sizeType c=off[r]; c<off[r+1]; c++)
                coeff+=SQR(row[c].first);
            coeff=std::sqrt(coeff);
            if(std::abs(coeff) < thres)
                coeff=0.0f;
            else coeff=1.0f/coeff;
#undef SQR
        }
    }
    static void scale_row(const vector_type& s,sparse_matrix_type& x) {
        const std::vector<sizeType>& off=x.getRowOffset();
        typename sparse_matrix_type::ROW& row=x.getValue();
        for(sizeType r=0; r<x.rows(); r++)
            for(sizeType c=off[r]; c<off[r+1]; c++)
                row[c].first*=s[r];
    }
    //matrix builder
    static void clear(matrix_builder& builder) {
        builder.clear();
    }
    static void add(size_t r,size_t c,value_type a,matrix_builder& builder) {
        builder.push_back(Eigen::Triplet<value_type,sizeType>(r,c,a));
    }
    static void set(size_t r,size_t c,matrix_builder& builder,sparse_matrix_type& x) {
        if(!x.getValue().empty()) {
            //in place set
            std::sort(builder.begin(),builder.end(),sparse_matrix_type::Lss());
            sizeType rows=x.rows(),rowB,rowE,c;
            typename matrix_builder::const_iterator itrB=builder.begin(),itrE=builder.end();
            for(sizeType r=0; r<rows; r++) {
                for(rowB=x.getRowOffset()[r],rowE=x.getRowOffset()[r+1]; rowB<rowE; rowB++) {
                    c=x.getValue()[rowB].second;
                    value_type& val=x.getValue()[rowB].first;
                    val=0.0f;
                    if(itrB < itrE && (itrB->row() != r || itrB->col() != c))
                        throw "in place set failure";
                    while(itrB < itrE && itrB->row() == r && itrB->col() == c) {
                        val+=itrB->value();
                        itrB++;
                    }
                }
            }
        } else {
            //reallocate and set
            x.resize(r,c);
            x.buildFromTripletsDepulicate(builder);
        }
    }
};

}
}

#endif
