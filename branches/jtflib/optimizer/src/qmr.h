#ifndef JTF_QMR_H
#define JTF_QMR_H

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>

#include <mpi.h>
#include <petscksp.h>
#include <petscversion.h>
#include <zjucad/linear_solver/dl_petsc/petsc_ver.h>

namespace jtf{

  // dense version
  class qmr_solver{
  public:
    ///
    /// @brief qmr_solver solve Ax=b,
    /// @param A
    /// @param b
    ///
    qmr_solver(const zjucad::matrix::matrix<double> & A,
               const zjucad::matrix::matrix<double> & b)
      :A_(A), b_(b), eps_(1e-7){}

    ///
    /// @brief The Quasi Minimal Residual Method without Look-ahead
    /// @param x
    /// @param max_iter
    /// @return       0 if work fine, or non-zeros
    ///
    int solve(double *x, size_t max_iter)const;

  private:
    const zjucad::matrix::matrix<double> &A_;
    const zjucad::matrix::matrix<double> &b_;
    const double eps_;
  };

  // dense version
  class qmr_solver_sparse{
  public:
    ///
    /// @brief qmr_solver solve Ax=b,
    /// @param A
    /// @param b
    ///
    qmr_solver_sparse(const hj::sparse::csc<double,int32_t> & A,
                      const zjucad::matrix::matrix<double> & b)
      :A_(A), b_(b), eps_(1e-6){}

    ///
    /// @brief The Quasi Minimal Residual Method without Look-ahead
    /// @param x
    /// @param max_iter
    /// @return       0 if work fine, or non-zeros
    ///
    int solve(double *x, size_t max_iter)const;

  private:
    const hj::sparse::csc<double,int32_t> & A_;
    const zjucad::matrix::matrix<double> &b_;
    const double eps_;
  };

  class PETsc_init
  {
  public:
    PETsc_init() {
      if(!PetscInitializeCalled){
          PetscInitializeNoArguments();
        }
    }
    virtual ~PETsc_init(){
      if(!PetscFinalizeCalled)
        PetscFinalize();
    }
  };

  // PETSC version TFQMR
  class PETsc_TFQMR {
  public:
    PETsc_TFQMR(bool transA,
                const hj::sparse::csc<double,int32_t> & A,
                const zjucad::matrix::matrix<double> & b)
      :eps_(1e-7)
    {
      comm = MPI_COMM_SELF;

      if(!transA){
          AT_ptr.reset(new hj::sparse::csc<double,int32_t>(A));
          hj::sparse::trans(A, *AT_ptr);
        }else{
          AT_ptr.reset(new hj::sparse::csc<double,int32_t>(A));
        }

      Dim = AT_ptr->size(1);
      if(sizeof(int32_t) == sizeof(int)){
          MatCreateSeqAIJWithArrays(
                comm, AT_ptr->size(2), AT_ptr->size(1),
                const_cast<int*> (&(AT_ptr->ptr()[0])),
              const_cast<int*> (&(AT_ptr->idx()[0])),
              &(AT_ptr->val()[0]), &A_);
        }else{
          ptr_ = AT_ptr->ptr();
          idx_ = AT_ptr->idx();
          MatCreateSeqAIJWithArrays(
                comm, AT_ptr->size(2), AT_ptr->size(1),
                &ptr_[0], &idx_[0], &(AT_ptr->val()[0]), &A_);
        }

      KSPCreate(comm,&solver);
      KSPSetOperators(solver, A_, A_, SAME_PRECONDITIONER);

      KSPGetPC(solver,&pc);
      PCSetType(pc,PCNONE);
      KSPSetFromOptions(solver);
      KSPSetType(solver, KSPMINRES);

      KSPSetInitialGuessNonzero(solver,PETSC_TRUE);
      KSPSetUp(solver);
      vec_create_seq_with_array(comm, 1, Dim, &b[0], &b_);
    }

    int solve(double *x, size_t max_iter){
      KSPSetTolerances(solver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
      vec_create_seq_with_array(comm, 1, Dim, x, &X_);
      PetscErrorCode err = KSPSolve(solver, b_, X_);

      vec_destroy(&b_);
      vec_destroy(&X_);
      return 0;
    }
    ~PETsc_TFQMR(){
      mat_destroy(&A_);
      KSP_destroy(&solver);
      //  PetscFinalize();
    }
  private:
    std::shared_ptr<hj::sparse::csc<double,int32_t> > AT_ptr;
    Mat A_;
    Vec b_, X_;
    zjucad::matrix::matrix<int> ptr_, idx_;
    MPI_Comm comm;
    KSP solver;
    PC pc;
    int Dim;
    const double eps_;
  };
}
#endif
