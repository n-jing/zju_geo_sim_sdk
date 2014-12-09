#ifndef _KERNEL_ZJUCAD_MATRIX_H_
#define _KERNEL_ZJUCAD_MATRIX_H_

#include <stdint.h>

#include <zjucad/matrix/iterator.h>
#include <zjucad/matrix/configure.h>
#include <zjucad/matrix/colon.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/operation.h>

#include <hjlib/sparse/format.h>
#include <hjlib/sparse/operation.h>

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename SCALAR_TYPE,typename INT_TYPE>
struct kernel_traits<hj::sparse::csc<SCALAR_TYPE,INT_TYPE> > {
    //typedefs
    typedef SCALAR_TYPE value_type;
    typedef INT_TYPE int_type;
    typedef typename matrix::matrix<value_type> vector_type;
    typedef typename hj::sparse::csc<value_type,INT_TYPE> sparse_matrix_type;
    //matrix builder
    typedef std::pair<INT_TYPE,INT_TYPE> coord;
    typedef std::pair<coord,value_type> triplet;
    typedef std::vector<triplet> matrix_builder;
    struct matrix_builder_less {
        bool operator()(const triplet& a,const triplet& b) const {
            return a.first.second < b.first.second || (a.first.second == b.first.second && a.first.first < b.first.first);
        }
    };
    //functions of blas: zjucad matrix implementation
    static value_type dot(const vector_type& x,const vector_type& y,int64_t n=-1) {
        if(n == -1)return matrix::dot(x,y);
        else return matrix::dot(x(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void()),
                                    y(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void()));
    }
    static void axpy(value_type a,const vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1)y+=x*a;
        else y(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void())+=
                x(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void())*a;
    }
    static void copy(const vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1)y=x;
        else y(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void())=
                x(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void());
    }
    static void swap(vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1) {
            assert(x.size(1) == y.size(1));
            n=x.size(1);
        }
        for(size_t i=0; i<(size_t)n; i++)
            std::swap(x[i],y[i]);
    }
    static value_type nrm2(const vector_type& x,int64_t n=-1) {
        if(n == -1)return matrix::norm(x);
        else return matrix::norm(x(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void()));
    }
    static value_type asum(const vector_type& x,int64_t n=-1) {
        if(n == -1)n=x.size(1);
        typename vector_type::const_iterator i=x.begin(), end=x.begin()+n;
        value_type r=0;
        for(; i!=end; ++i)r+=std::abs(*i);
        return r;
    }
    static void scal(value_type a,vector_type& x,int64_t n=-1) {
        if(n == -1)x*=a;
        else x(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void())*=a;
    }
    static value_type amax(const vector_type& x,int64_t n=-1) {
        if(n == -1)n=x.size(1);
        typename vector_type::const_iterator i=x.begin(), end=x.begin()+n;
        value_type r=0;
        for(; i!=end; ++i)r=std::max<value_type>(r,std::abs(*i));
        return r;
    }
    //other necessary functions
    static void resize(size_t n,vector_type& x) {
        x.resize(n,1);
    }
    static void zero(vector_type& x,int64_t n=-1) {
        set(0.0f,x,n);
    }
    static void set(value_type a,vector_type& x,int64_t n=-1) {
        if(n == -1)x=matrix::ones<value_type>(x.size(1),1)*a;
        else x(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void())=matrix::ones<value_type>(n,1)*a;
    }
    static void set(size_t i,value_type a,vector_type& x) {
        x[i]=a;
    }
    static value_type get(size_t i,const vector_type& x) {
        return x[i];
    }
    static size_t rows(const vector_type& x) {
        return x.size(1);
    }
    static void cmul(const vector_type& x,vector_type& y,int64_t n=-1) {
        if(n == -1)y=ele_prod(x,y);
        else y(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void())=
                ele_prod(x(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void()),
                         y(matrix::colon<int64_t,int64_t>(0,n-1),matrix::colon_void()));
    }
    //matrix vector multiplication function
    static void mul(const sparse_matrix_type& A,const vector_type& x,vector_type& Ax) {
        zero(Ax,rows(A));
        mul_add(A,x,Ax);
    }
    static void mul_add(const sparse_matrix_type& A,const vector_type& x,vector_type& bPAx) {
        hj::sparse::mv(false,A,x,bPAx);
    }
    static void mul_sub(const sparse_matrix_type& A,const vector_type& x,vector_type& bSAx) {
        vector_type tmp;
        resize(rows(A),tmp);
        mul(A,x,tmp);
        axpy(-1.0f,tmp,bSAx,rows(A));
    }
    static void mul_t(const sparse_matrix_type& A,const vector_type& x,vector_type& ATx) {
        zero(ATx,cols(A));
        mul_t_add(A,x,ATx);
    }
    static void mul_t_add(const sparse_matrix_type& A,const vector_type& x,vector_type& bPATx) {
        hj::sparse::mv(true,A,x,bPATx);
    }
    static void mul_t_sub(const sparse_matrix_type& A,const vector_type& x,vector_type& bSATx) {
        vector_type tmp;
        resize(cols(A),tmp);
        mul_t(A,x,tmp);
        axpy(-1.0f,tmp,bSATx,cols(A));
    }
    //other necessary matrix functions
    static void resize(size_t r,size_t c,sparse_matrix_type& x) {
        x.resize((INT_TYPE)r,(INT_TYPE)c,0);
    }
    static void set(size_t r,size_t c,value_type a,sparse_matrix_type& x) {
        //check for in-place put
        INT_TYPE i=x.ptr()[c],j;
        for(; i<x.ptr()[c+1]; i++) {
            if(x.idx()[i] == r) {
                x.val()[i]=a;
                return;
            } else if(x.idx()[i] > (INT_TYPE)r)
                break;
        }

        //set ptr_
        for(j=c+1; j<x.ptr().size(1); j++)
            x.ptr()[j]++;

        //set idx_ and val_ after c
        matrix::matrix<INT_TYPE> idx_new;
        idx_new.resize(x.idx().size(1)+1,1);
        if(i > 0)
            idx_new(matrix::colon<int64_t,int64_t>(0,i-1),matrix::colon_void())=
                x.idx()(matrix::colon<int64_t,int64_t>(0,i-1),matrix::colon_void());
        if(i < x.idx().size(1))
            idx_new(matrix::colon<int64_t,int64_t>(i+1,idx_new.size(1)-1),matrix::colon_void())=
                x.idx()(matrix::colon<int64_t,int64_t>(i,x.idx().size(1)-1),matrix::colon_void());
        idx_new[i]=r;

        matrix::matrix<value_type> val_new;
        val_new.resize(x.val().size(1)+1,1);
        if(i > 0)
            val_new(matrix::colon<int64_t,int64_t>(0,i-1),matrix::colon_void())=
                x.val()(matrix::colon<int64_t,int64_t>(0,i-1),matrix::colon_void());
        if(i < x.val().size(1))
            val_new(matrix::colon<int64_t,int64_t>(i+1,val_new.size(1)-1),matrix::colon_void())=
                x.val()(matrix::colon<int64_t,int64_t>(i,x.val().size(1)-1),matrix::colon_void());
        val_new[i]=a;

        x.idx()=idx_new;
        x.val()=val_new;
    }
    static size_t rows(const sparse_matrix_type& x) {
        return x.size(1);
    }
    static size_t cols(const sparse_matrix_type& x) {
        return x.size(2);
    }
    static void row_scale_coeff(const sparse_matrix_type& x,vector_type& s,value_type thres) {
        zero(s,rows(x));
#define SQR(x) (x)*(x)
        for(INT_TYPE c=0; c<(INT_TYPE)x.val().size(); c++)
            s[x.idx()[c]]+=SQR(x.val()[c]);
        for(INT_TYPE c=0; c<(INT_TYPE)x.size(1); c++) {
            s[c]=std::sqrt(s[c]);
            if(std::abs(s[c]) < thres)
                s[c]=0.0f;
            else s[c]=1.0f/s[c];
        }
#undef SQR
    }
    static void scale_row(const vector_type& s,sparse_matrix_type& x) {
        for(INT_TYPE c=0; c<x.val().size(); c++)
            x.val()[c]*=s[x.idx()[c]];
    }
    //matrix builder
    static void clear(matrix_builder& builder) {
        builder.clear();
    }
    static void add(size_t r,size_t c,value_type a,matrix_builder& builder) {
        builder.push_back(triplet(coord((INT_TYPE)r,(INT_TYPE)c),a));
    }
    static void set(size_t r,size_t c,matrix_builder& builder,sparse_matrix_type& x) {
        std::sort(builder.begin(),builder.end(),matrix_builder_less());
        if(x.idx().size() > 0) {
            //in place set
            INT_TYPE cols=(INT_TYPE)x.size(2),colB,colE;
            typename matrix_builder::const_iterator itrB=builder.begin(),itrE=builder.end();
            for(coord c(0,0); c.second<cols; c.second++) {
                for(colB=x.ptr()[c.second],colE=x.ptr()[c.second+1]; colB<colE; colB++) {
                    c.first=x.idx()[colB];
                    value_type& val=x.val()[colB];
                    val=0.0f;
                    if(itrB < itrE && itrB->first != c)
                        throw "in place set failure";
                    while(itrB < itrE && itrB->first == c) {
                        val+=itrB->second;
                        itrB++;
                    }
                }
            }
        } else {
            //reallocate and set
            resize(r,c,x);
            size_t i=0,n=builder.size(),nn=n-1,lc=0;
            x.ptr()[0]=0;
            std::vector<INT_TYPE> idx;
            std::vector<value_type> vals;
            idx.reserve(builder.size());
            vals.reserve(builder.size());
            for(; i<n; i++) {
                while(lc < (size_t)(builder[i].first.second)) {
                    lc++;
                    x.ptr()[lc]=(INT_TYPE)idx.size();
                }
                value_type val=builder[i].second;
                while(i < nn && builder[i].first == builder[i+1].first) {
                    i++;
                    val+=builder[i].second;
                }
                idx.push_back(builder[i].first.first);
                vals.push_back(val);
            }
            while(lc < c) {
                lc++;
                x.ptr()[lc]=(INT_TYPE)idx.size();
            }
            x.idx().resize((INT_TYPE)idx.size());
            x.val().resize((INT_TYPE)idx.size());
            std::copy(idx.begin(),idx.end(),x.idx().begin());
            std::copy(vals.begin(),vals.end(),x.val().begin());
        }
    }
};

}
}

#endif
