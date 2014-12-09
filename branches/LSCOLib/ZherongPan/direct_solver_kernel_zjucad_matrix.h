#ifndef _DIRECT_SOLVER_KERNEL_ZJUCAD_MATRIX_H_
#define _DIRECT_SOLVER_KERNEL_ZJUCAD_MATRIX_H_

#include <boost/property_tree/ptree.hpp>
#include <hjlib/sparse/operation.h>
#include <hjlib/sparse/format.h>
#include <Eigen/SPQRSupport>
#include <Eigen/SparseQR>
#include <Eigen/Sparse>

namespace zjucad
{
namespace LSCO
{

template <typename MATRIX_TYPE>
struct kernel_traits;

template <typename MATRIX_TYPE>
struct direct_solver_kernel;

template <typename SCALAR_TYPE,typename INT_TYPE>
class direct_solver_kernel<typename hj::sparse::csc<SCALAR_TYPE,INT_TYPE> >
{
public:
    //typedefs
    typedef kernel_traits<typename hj::sparse::csc<SCALAR_TYPE,INT_TYPE> > kernel_type;
    typedef typename kernel_type::value_type value_type;
    typedef typename kernel_type::vector_type vector_type;
    typedef typename kernel_type::sparse_matrix_type sparse_matrix_type;
    //functions
    direct_solver_kernel(boost::property_tree::ptree& opt):m_opt(opt) {
        m_opt.add_child("direct_solver_in",boost::property_tree::ptree());
        m_opt.add_child("direct_solver_out",boost::property_tree::ptree());
        m_opt.put<int64_t>("direct_solver_out.n_pos",-1);
        m_opt.put<int64_t>("direct_solver_out.n_neg",-1);
        m_opt.put<int64_t>("direct_solver_out.n_zero",-1);
        m_sol.setPivotThreshold(m_opt.get<value_type>("direct_solver_in.threshold",1E-6f));
    }
    size_t rows() const {
        if(m_use_faster)
            return (size_t)(m_sol_faster.rows());
        else
            return (size_t)(m_sol.rows());
    }
    void build_AAT(const sparse_matrix_type& A,bool AAT,bool non_singular ) {
        Eigen::MappedSparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> m
        ((INT_TYPE)A.size(1),
         (INT_TYPE)A.size(2),
         (INT_TYPE)hj::sparse::nnz(A),
         (INT_TYPE*) &(A.ptr_[0]),
         (INT_TYPE*)&(A.idx_[0]),
         const_cast<SCALAR_TYPE*>(&(A.val_[0])));
        Eigen::SparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> P_mat;
        if(AAT) P_mat = m*m.transpose();
        else P_mat = m.transpose() * m;

        if(non_singular) {
            m_use_faster = true;
            m_sol_faster.compute(P_mat);
        } else {
            m_use_faster = false;
            m_sol.compute(P_mat);
        }
        get_info();
    }
    void build_H(const sparse_matrix_type& H,bool non_singular,bool positive) {
        Eigen::SparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> P_mat;
        to_eigen(H, P_mat);

        if(positive) {
            m_use_faster = true;
            m_sol_faster.compute(P_mat);
        } else {
            m_use_faster = false;
            m_sol.compute(P_mat);
        }
        get_info();
    }
    void build_KKT(const sparse_matrix_type& H,const sparse_matrix_type& C,bool non_singular) {
        std::vector<Eigen::Triplet<SCALAR_TYPE,INT_TYPE> > trips;
        to_eigen(H,0,0,trips);
        to_eigen(C,(INT_TYPE)H.size(1),0,trips);
        to_eigen(C,0,(INT_TYPE)H.size(2),trips, true);

        Eigen::SparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> KKT;
        KKT.resize((INT_TYPE)(H.size(1)+C.size(1)), (INT_TYPE)(H.size(1)+C.size(1)));
        KKT.setFromTriplets(trips.begin(), trips.end());

        m_use_faster = false;
        m_sol.compute(KKT);
        get_info();
    }
    void build_KKT_schur(const sparse_matrix_type& invH,const sparse_matrix_type& C,bool non_singular) {
        Eigen::SparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> invH_mat,C_mat;
        to_eigen(invH, invH_mat);
        to_eigen(C, C_mat);
        if(non_singular) {
            m_use_faster = true;
            m_sol_faster.compute(C_mat*invH_mat*C_mat.transpose());
        } else {
            m_use_faster = false;
            m_sol.compute(C_mat*invH_mat*C_mat.transpose());
        }
        get_info();
    }
    void solve(const vector_type& rhs,vector_type& result) {
        Eigen::Map<Eigen::Matrix<value_type,-1,-1> > rx(&result[0], result.size(1),result.size(2));
        Eigen::Map<const Eigen::Matrix<value_type,-1,-1> > rb(&rhs[0], rhs.size(1), rhs.size(2));
        //m_use_faster = false;
        if(m_use_faster)
            rx = m_sol_faster.solve(rb);
        else {
            rx = m_sol.solve(rb);
        }
    }
private:
    void get_info() {
        if(!m_use_faster) {
            if(m_sol.info() == Eigen::Success) {
                m_opt.put<std::string>("direct_solver_out.termination_type","successful");
                m_opt.put<int64_t>("direct_solver_out.n_zero",rows()-m_sol.rank());
            } else {
                m_opt.put("direct_solver_out.termination_type","fail");
                m_opt.put("direct_solver_out.n_zero",-1);
            }
        } else {
            if(m_sol_faster.info() == Eigen::Success) {
                m_opt.put<std::string>("direct_solver_out.termination_type","successful");
            } else {
                m_opt.put<std::string>("direct_solver_out.termination_type","fail");
            }
        }
    }
    void to_eigen(const sparse_matrix_type &A, Eigen::SparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> &a) {
        Eigen::MappedSparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> m
        ((INT_TYPE)A.size(1),
         (INT_TYPE)A.size(2),
         (INT_TYPE)hj::sparse::nnz(A),
         (INT_TYPE*) &(A.ptr_[0]),
         (INT_TYPE*)&(A.idx_[0]),
         const_cast<SCALAR_TYPE*>(&(A.val_[0])));
        a = m;
    }
    void to_eigen(const sparse_matrix_type &A, INT_TYPE pre_rows, INT_TYPE pre_cols,std::vector<Eigen::Triplet<SCALAR_TYPE,INT_TYPE> > & trips,bool trans = false) {
        if(trans == false) {
            for(INT_TYPE i = 0; i < A.ptr_.size()-1; ++i) {
                for(INT_TYPE j = A.ptr_[i]; j < A.ptr_[i+1]; ++j) {
                    trips.push_back(Eigen::Triplet<SCALAR_TYPE,INT_TYPE>(A.idx_[j]+pre_rows,i+pre_cols,A.val_[j]));
                }
            }
        } else {
            for(INT_TYPE i = 0; i < A.ptr_.size()-1; ++i) {
                for(INT_TYPE j = A.ptr_[i]; j < A.ptr_[i+1]; ++j) {
                    trips.push_back(Eigen::Triplet<SCALAR_TYPE,INT_TYPE>(i+pre_rows, A.idx_[j]+pre_cols,A.val_[j]));
                }
            }
        }
    }
private:
    boost::property_tree::ptree& m_opt;
    Eigen::SparseQR<Eigen::SparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> , Eigen::AMDOrdering<INT_TYPE> > m_sol;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<SCALAR_TYPE,Eigen::ColMajor,INT_TYPE> > m_sol_faster;
    bool m_use_faster;
};

}
}

#endif
