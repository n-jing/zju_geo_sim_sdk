#ifndef _DIRECT_SOLVER_KERNEL_H_
#define _DIRECT_SOLVER_KERNEL_H_

#include <stdint.h>
#include <cmath>
#include <boost/shared_array.hpp>
#include <boost/property_tree/ptree.hpp>

extern "C" {
    void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], double d[], int *it_num, int *rot_num );
}

namespace zjucad
{
namespace LSCO
{

//accepted parameter,default value
//direct_solver_in.dynamic,true
//direct_solver_in.require_inertia,false
//direct_solver_in.threshold,1E-6f

//returned parameter
//direct_solver_out.n_pos(-1 for unknown)
//direct_solver_out.n_neg(-1 for unknown)	//interior/direct require the factorization to be rank revealing
//direct_solver_out.n_zero(-1 for unknown)
//direct_solver_out.termination_type,std::string
//	successful
//	fail

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename MATRIX_TYPE>
struct direct_solver_kernel;

template <typename SCALAR_TYPE>
struct direct_solver_kernel<plain_matrix<SCALAR_TYPE> > {
public:
    //typedefs
    typedef kernel_traits<plain_matrix<SCALAR_TYPE> > kernel_type;
    typedef typename kernel_type::value_type value_type;
    typedef typename kernel_type::vector_type vector_type;
    typedef typename kernel_type::sparse_matrix_type sparse_matrix_type;
    //functions
    direct_solver_kernel(boost::property_tree::ptree& opt):m_opt(opt) {
        m_opt.add_child("direct_solver_in",boost::property_tree::ptree());
        m_opt.add_child("direct_solver_out",boost::property_tree::ptree());
    }
    size_t rows() const {
        return kernel_type::rows(m_ev);
    }
    ///
    /// \brief build  Px=b
    /// \param A      input csc matrix A
    /// \param AAT    if AAT=true, P=AAT, or P=ATA
    ///
    void build_AAT(const sparse_matrix_type& A,bool AAT,bool non_singular) {
        sparse_matrix_type AA;
        if(AAT)abt(A,A,AA);
        else atb(A,A,AA);
        build(AA);
    }
    ///
    /// \brief build        Hx=b
    /// \param H            input csc matrix H
    /// \param non_singular  whether H is nosingular or not
    /// \param positive     whether H is positive semi-definite or not
    ///
    void build_H(const sparse_matrix_type& H,bool non_singular,bool positive) {
        build(H);
    }
    ///
    /// \brief build        H CT  x =  b
    ///                     C 0
    /// \param H            input csc matrix H
    /// \param C            input csc matrix C
    /// \param non_singular  whether (H CT) is nonsingular or not
    ///                              (C  0)
    ///
    void build_KKT(const sparse_matrix_type& H,const sparse_matrix_type& C,bool non_singular) {
        assert(kernel_type::cols(C) == kernel_type::rows(H));
        assert(kernel_type::rows(H) == kernel_type::cols(H));

        sparse_matrix_type KKTm;
        kernel_type::resize(kernel_type::rows(H)+kernel_type::rows(C),
                            kernel_type::rows(H)+kernel_type::rows(C),KKTm);

        //H
        for(size_t r=0; r<kernel_type::rows(H); r++)
            for(size_t c=0; c<kernel_type::rows(H); c++)
                KKTm.m_data[r][c]=H.m_data[r][c];

        //C
        for(size_t r=0; r<kernel_type::rows(C); r++)
            for(size_t c=0; c<kernel_type::rows(H); c++) {
                KKTm.m_data[r+kernel_type::rows(H)][c]=C.m_data[r][c];
                KKTm.m_data[c][r+kernel_type::rows(H)]=C.m_data[r][c];
            }

        //zero
        for(size_t r=0; r<kernel_type::rows(C); r++)
            for(size_t c=0; c<kernel_type::rows(C); c++)
                KKTm.m_data[kernel_type::rows(H)+r][kernel_type::rows(H)+c]=0.0f;

        build(KKTm);
    }
    ///
    /// \brief build        (C*invH*C^T) x=b
    /// \param H            input csc matrix invH and csc matrix C
    /// \param non_singular whether C*invH*C^T is nosingular or not
    ///
    void build_KKT_schur(const sparse_matrix_type& invH,const sparse_matrix_type& C,bool non_singular) {
        sparse_matrix_type invH_CT,C_invH_CT;
        abt(invH,C,invH_CT);
        ab(C,invH_CT,C_invH_CT);
        build(C_invH_CT);
    }
    ///
    /// \brief solve      solve Px=b
    /// \param rhs        input b
    /// \param result     output result
    ///
    void solve(const vector_type& rhs,vector_type& result) {
        vector_type tmp;
        kernel_type::resize(rows(),tmp);
        kernel_type::mul_t(m_ev,rhs,tmp);
        kernel_type::cmul(m_e,tmp);
        kernel_type::resize(kernel_type::rows(rhs),result);
        kernel_type::mul(m_ev,tmp,result);
    }
private:
    static void abt(const sparse_matrix_type& a,const sparse_matrix_type& b,sparse_matrix_type& abt) {
        assert(kernel_type::cols(a) == kernel_type::cols(b));
        kernel_type::resize(kernel_type::rows(a),kernel_type::rows(b),abt);
        for(size_t r=0; r<kernel_type::rows(a); r++)
            for(size_t c=0; c<kernel_type::rows(b); c++) {
                abt.m_data[r][c]=0.0f;
                for(size_t d=0; d<kernel_type::cols(a); d++)
                    abt.m_data[r][c]+=a.m_data[r][d]*b.m_data[c][d];
            }
    }
    static void atb(const sparse_matrix_type& a,const sparse_matrix_type& b,sparse_matrix_type& atb) {
        assert(kernel_type::rows(a) == kernel_type::rows(b));
        kernel_type::resize(kernel_type::cols(a),kernel_type::cols(b),atb);
        for(size_t r=0; r<kernel_type::cols(a); r++)
            for(size_t c=0; c<kernel_type::cols(b); c++) {
                atb.m_data[r][c]=0.0f;
                for(size_t d=0; d<kernel_type::rows(a); d++)
                    atb.m_data[r][c]+=a.m_data[d][r]*b.m_data[d][c];
            }
    }
    static void ab(const sparse_matrix_type& a,const sparse_matrix_type& b,sparse_matrix_type& ab) {
        assert(kernel_type::cols(a) == kernel_type::rows(b));
        kernel_type::resize(kernel_type::rows(a),kernel_type::cols(b),ab);
        for(size_t r=0; r<kernel_type::rows(a); r++)
            for(size_t c=0; c<kernel_type::cols(b); c++) {
                ab.m_data[r][c]=0.0f;
                for(size_t d=0; d<kernel_type::cols(a); d++)
                    ab.m_data[r][c]+=a.m_data[r][d]*b.m_data[d][c];
            }
    }
    void build(const sparse_matrix_type& m) {
        m_opt.get_child("direct_solver_out").clear();

        int n=(int)kernel_type::rows(m);
        boost::shared_array<double> a(new double[n*n]);
        for(int64_t r=0; r<n; r++)
            for(int64_t c=0; c<n; c++)
                a[r*n+c]=m.m_data[r][c];
        boost::shared_array<double> d(new double[n]);
        boost::shared_array<double> v(new double[n*n]);

        int it_num,rot_num;
        jacobi_eigenvalue(n,a.get(),10000,v.get(),d.get(),&it_num,&rot_num);

        kernel_type::resize(n,m_e);
        kernel_type::resize(n,n,m_ev);
        int n_pos=0,n_neg=0,n_zero=0;
        for(int i=0; i<n; i++) {
            if(std::abs(d[i]) < m_opt.get<value_type>("direct_solver_in.threshold",1E-6f)) {
                m_e.m_data[i]=0.0f;
                n_zero++;
            } else if(d[i] > 0.0f) {
                m_e.m_data[i]=1.0f/(value_type)d[i];
                n_pos++;
            } else {
                m_e.m_data[i]=1.0f/(value_type)d[i];
                n_neg++;
            }
            for(int j=0; j<n; j++)
                m_ev.m_data[j][i]=(value_type)v[i*n+j];
        }

        m_opt.put<int64_t>("direct_solver_out.n_pos",n_pos);
        m_opt.put<int64_t>("direct_solver_out.n_neg",n_neg);
        m_opt.put<int64_t>("direct_solver_out.n_zero",n_zero);
        m_opt.put<std::string>("direct_solver_out.termination_type","successful");
    }
    boost::property_tree::ptree& m_opt;
    sparse_matrix_type m_ev;
    vector_type m_e;
};

}
}

#endif
