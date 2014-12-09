#ifndef _DIRECT_SOLVER_KERNEL_COMMON_FILE_PZR_H_
#define _DIRECT_SOLVER_KERNEL_COMMON_FILE_PZR_H_

#include <stdint.h>
#include <boost/property_tree/ptree.hpp>
#include <SuiteSparseQR.hpp>
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

template <typename SCALAR_TYPE>
class direct_solver_kernel<COMMON::FixedSparseMatrix<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > >
{
public:
    //typedefs
    typedef kernel_traits<COMMON::FixedSparseMatrix<SCALAR_TYPE,COMMON::Kernel<SCALAR_TYPE> > > kernel_type;
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
    void build_AAT(const sparse_matrix_type& A,bool AAT,bool non_singular) {
        Eigen::SparseMatrix<value_type,0,sizeType> AE,HE;
        A.toEigen(AE);
        if(AAT)
            HE=AE*AE.transpose();
        else
            HE=AE.transpose()*AE;
        if(non_singular) {
            m_use_faster=true;
            m_sol_faster.compute(HE);
        } else {
            m_use_faster=false;
            m_sol.compute(HE);
        }
        get_info();
    }
    void build_H(const sparse_matrix_type& H,bool non_singular,bool positive) {
        Eigen::SparseMatrix<value_type,0,sizeType> HE;
        H.toEigen(HE);
        if(positive) {
            m_use_faster=true;
            m_sol_faster.compute(HE);
        } else {
            m_use_faster=false;
            m_sol.compute(HE);
        }
        get_info();
    }
    void build_KKT(const sparse_matrix_type& H,const sparse_matrix_type& C,bool non_singular) {
        std::vector<Eigen::Triplet<value_type,sizeType> > trips;
        H.toEigen(0,0,trips);
        C.toEigen(H.rows(),0,trips,true);

        Eigen::SparseMatrix<value_type,0,sizeType> KKTE;
        KKTE.resize((int)(H.rows()+C.rows()),(int)(H.rows()+C.rows()));
        KKTE.setFromTriplets(trips.begin(),trips.end());

        m_use_faster=false;
        m_sol.compute(KKTE);
        get_info();
    }
    void build_KKT_schur(const sparse_matrix_type& invH,const sparse_matrix_type& C,bool non_singular) {
        Eigen::SparseMatrix<value_type,0,sizeType> invH_mat,C_mat;
        invH.toEigen(invH_mat);
        C.toEigen(C_mat);

        if(non_singular) {
            m_use_faster=true;
            m_sol_faster.compute(C_mat*invH_mat*C_mat.transpose());
        } else {
            m_use_faster=false;
            m_sol.compute(C_mat*invH_mat*C_mat.transpose());
        }
        get_info();
    }
    void solve(const vector_type& rhs,vector_type& result) {
        if(m_use_faster)
            result=m_sol_faster.solve(rhs);
        else result=m_sol.solve(rhs);
    }
private:
    void get_info() {
        if(!m_use_faster) {
            if(m_sol.info() == Eigen::Success) {
                m_opt.put<std::string>("direct_solver_out.termination_type","successful");
                m_opt.put<int64_t>("direct_solver_out.n_zero",rows()-m_sol.rank());
            } else {
                m_opt.put<std::string>("direct_solver_out.termination_type","fail");
                m_opt.put<int64_t>("direct_solver_out.n_zero",-1);
            }
        } else {
            if(m_sol_faster.info() == Eigen::Success) {
                m_opt.put<std::string>("direct_solver_out.termination_type","successful");
            } else {
                m_opt.put<std::string>("direct_solver_out.termination_type","fail");
            }
        }
    }
    boost::property_tree::ptree& m_opt;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<value_type,0,sizeType> > m_sol_faster;
    Eigen::SparseQR<Eigen::SparseMatrix<value_type,0,sizeType>,Eigen::COLAMDOrdering<sizeType> > m_sol;
    bool m_use_faster;
};

}
}

#endif