#ifndef _KERNEL_H_
#define _KERNEL_H_

#include <stdint.h>
#include <vector>
#include <cstddef>
#include <cmath>
#include <assert.h>

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename SCALAR_TYPE=double>
struct plain_vector {
    std::vector<SCALAR_TYPE> m_data;
};
template <typename SCALAR_TYPE=double>
struct plain_matrix {
    //row then col
    std::vector<std::vector<SCALAR_TYPE> > m_data;
};
template <typename SCALAR_TYPE>
struct kernel_traits<plain_matrix<SCALAR_TYPE> > {
    //typedefs
    typedef SCALAR_TYPE value_type;
    typedef plain_vector<value_type> vector_type;
    typedef plain_matrix<value_type> sparse_matrix_type;
    struct matrix_builder {
        std::vector<size_t> m_row;
        std::vector<size_t> m_col;
        std::vector<value_type> m_val;
    };
    //functions of blas: plain slow implementation
    static value_type dot(const vector_type& x,const vector_type& y,int64_t n=-1) {
        if(n == -1) {
            assert(x.m_data.size() == y.m_data.size());
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n && y.m_data.size() >= (size_t)n);
        }
        value_type ret=0.0f;
        for(size_t i=0; i<(size_t)n; i++)
            ret+=x.m_data[i]*y.m_data[i];
        return ret;
    }
    static void axpy(value_type a,const vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1) {
            assert(x.m_data.size() == y.m_data.size());
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n && y.m_data.size() >= (size_t)n);
        }
        for(size_t i=0; i<(size_t)n; i++)
            y.m_data[i]+=x.m_data[i]*a;
    }
    static void copy(const vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1) {
            assert(x.m_data.size() == y.m_data.size());
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n && y.m_data.size() >= (size_t)n);
        }
        for(size_t i=0; i<(size_t)n; i++)
            y.m_data[i]=x.m_data[i];
    }
    static void swap(vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1) {
            assert(x.m_data.size() == y.m_data.size());
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n && y.m_data.size() >= (size_t)n);
        }
        for(size_t i=0; i<(size_t)n; i++)
            std::swap(y.m_data[i],x.m_data[i]);
    }
    static value_type nrm2(const vector_type& x,int64_t n=-1) {
        return sqrt(dot(x,x,n));
    }
    static value_type asum(const vector_type& x,int64_t n=-1) {
        if(n == -1) {
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n);
        }
        value_type ret=0.0f;
        for(size_t i=0; i<(size_t)n; i++)
            ret+=std::abs(x.m_data[i]);
        return ret;
    }
    static void scal(value_type a,vector_type& x,int64_t n=-1) {
        if(n == -1) {
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n);
        }
        for(size_t i=0; i<(size_t)n; i++)
            x.m_data[i]*=a;
    }
    static value_type amax(const vector_type& x,int64_t n=-1) {
        if(n == -1) {
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n);
        }
        value_type ret=0.0f;
        for(size_t i=0; i<(size_t)n; i++)
            ret=std::max<value_type>(ret,std::abs(x.m_data[i]));
        return ret;
    }
    //other necessary functions
    static void resize(size_t n,vector_type& x) {
        x.m_data.resize(n);
    }
    static void zero(vector_type& x,int64_t n=-1) {
        set(0.0f,x,n);
    }
    static void set(value_type a,vector_type& x,int64_t n=-1) {
        if(n == -1) {
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n);
        }
        for(size_t i=0; i<(size_t)n; i++)
            x.m_data[i]=a;
    }
    static void set(size_t i,value_type a,vector_type& x) {
        x.m_data[i]=a;
    }
    static value_type get(size_t i,const vector_type& x) {
        return x.m_data[i];
    }
    static size_t rows(const vector_type& x) {
        return x.m_data.size();
    }
    static void cmul(const vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1) {
            assert(x.m_data.size() == y.m_data.size());
            n=x.m_data.size();
        } else {
            assert(x.m_data.size() >= (size_t)n && y.m_data.size() >= (size_t)n);
        }
        for(size_t i=0; i<(size_t)n; i++)
            y.m_data[i]*=x.m_data[i];
    }
    //matrix vector multiplication function
    static void mul(const sparse_matrix_type& A,const vector_type& x,vector_type& Ax) {
        zero(Ax,rows(A));
        mul_add(A,x,Ax);
    }
    static void mul_add(const sparse_matrix_type& A,const vector_type& x,vector_type& bPAx) {
        for(size_t r=0; r<A.m_data.size(); r++)
            for(size_t c=0; c<A.m_data[r].size(); c++)
                bPAx.m_data[r]+=A.m_data[r][c]*x.m_data[c];
    }
    static void mul_sub(const sparse_matrix_type& A,const vector_type& x,vector_type& bSAx) {
        for(size_t r=0; r<A.m_data.size(); r++)
            for(size_t c=0; c<A.m_data[r].size(); c++)
                bSAx.m_data[r]-=A.m_data[r][c]*x.m_data[c];
    }
    static void mul_t(const sparse_matrix_type& A,const vector_type& x,vector_type& ATx) {
        zero(ATx,cols(A));
        mul_t_add(A,x,ATx);
    }
    static void mul_t_add(const sparse_matrix_type& A,const vector_type& x,vector_type& bPATx) {
        for(size_t r=0; r<A.m_data.size(); r++)
            for(size_t c=0; c<A.m_data[r].size(); c++)
                bPATx.m_data[c]+=A.m_data[r][c]*x.m_data[r];
    }
    static void mul_t_sub(const sparse_matrix_type& A,const vector_type& x,vector_type& bSATx) {
        for(size_t r=0; r<A.m_data.size(); r++)
            for(size_t c=0; c<A.m_data[r].size(); c++)
                bSATx.m_data[c]-=A.m_data[r][c]*x.m_data[r];
    }
    //other necessary matrix functions
    static void resize(size_t r,size_t c,sparse_matrix_type& x) {
        x.m_data.assign(r,std::vector<value_type>());
        for(size_t rr=0; rr<r; rr++)
            x.m_data[rr].assign(c,0.0f);
    }
    static void set(size_t r,size_t c,value_type a,sparse_matrix_type& x) {
        x.m_data[r][c]=a;
    }
    static size_t rows(const sparse_matrix_type& x) {
        return x.m_data.size();
    }
    static size_t cols(const sparse_matrix_type& x) {
        return x.m_data.empty() ? 0 : x.m_data[0].size();
    }
    static void row_scale_coeff(const sparse_matrix_type& x,vector_type& s,value_type thres) {
        for(size_t r=0; r<rows(x); r++) {
#define SQR(x) (x)*(x)
            value_type& coeff=s.m_data[r];
            coeff=0.0f;
            for(size_t c=0; c<cols(x); c++)
                coeff+=SQR(x.m_data[r][c]);
            coeff=std::sqrt(coeff);
            if(std::abs(coeff) < thres)
                coeff=0.0f;
            else coeff=1.0f/coeff;
#undef SQR
        }
    }
    static void scale_row(const vector_type& s,sparse_matrix_type& x) {
        for(size_t r=0; r<rows(x); r++)
            for(size_t c=0; c<cols(x); c++)
                x.m_data[r][c]*=s.m_data[r];
    }
    //matrix builder
    static void clear(matrix_builder& builder) {
        builder.m_row.clear();
        builder.m_col.clear();
        builder.m_val.clear();
    }
    static void add(size_t r,size_t c,value_type a,matrix_builder& builder) {
        builder.m_row.push_back(r);
        builder.m_col.push_back(c);
        builder.m_val.push_back(a);
    }
    static void set(size_t r,size_t c,matrix_builder& builder,sparse_matrix_type& x) {
        resize(r,c,x);
        for(size_t i=0; i<builder.m_row.size(); i++)
            x.m_data[builder.m_row[i]][builder.m_col[i]]+=builder.m_val[i];
    }
};

}
}

#endif
