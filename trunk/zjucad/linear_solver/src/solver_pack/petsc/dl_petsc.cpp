#include <iostream>
#include <petscksp.h>
#include "./include/petsc_ver.h"
#if PETSC_VERSION_GT(3, 2, 0)
#include <petsc-private/matimpl.h>
#else
#include <private/matimpl.h>
#endif
#include <petscversion.h>
#include "../../dl_petsc.h"
#include <hjlib/sparse/sparse.h>
#include <zjucad/matrix/itr_matrix.h>

class PETsc_imp : public PETsc
{
public:
  PETsc_imp() {
    std::cout << "call PETsc_imp" << std::endl;
    if(!PetscInitializeCalled){
        int argc = 0;
        //		PetscInitialize(&argc, 0,0,0);
        PetscInitializeNoArguments();
      }
  }
  virtual ~PETsc_imp() {
    if(!PetscFinalizeCalled){
        // std::cout << "PETsc finalized" << std::endl;
        PetscFinalize();
      }
  }
};// pets_init;

class PETsc_CG_imp : public PETsc_CG
{
public:
  PETsc_CG_imp(const double * val, const int32_t * idx,  const int32_t * ptr,
               const size_t nnz, const size_t row, const size_t col, const char *pc_str) {
    std::cout << "call PETsc_CG_imp" << std::endl;
    comm = MPI_COMM_SELF;
    Dim = row;
    if(sizeof(int32_t) == sizeof(int)) {
        MatCreateSeqAIJWithArrays(
              comm,
              row, col,
              const_cast<int *>(ptr), const_cast<int *>(idx), const_cast<double *>(val), &A_);
      }
    else {
        zjucad::matrix::itr_matrix<const int32_t*> ptr_m(col, 1, ptr);
        zjucad::matrix::itr_matrix<const int32_t*> idx_m(nnz, 1, idx);
        ptr_ = ptr_m;
        idx_ = idx_m;
        MatCreateSeqAIJWithArrays(
              comm,
              row, col, &ptr_[0], &idx_[0], const_cast<double *>(val), &A_);
      }
    KSPCreate(comm,&solver);
    KSPSetOperators(solver, A_, A_, SAME_PRECONDITIONER);
    PC pc;
    KSPGetPC(solver,&pc);
    PCSetType(pc, pc_str);
    KSPSetType(solver, KSPCG);
    KSPCGSetType(solver, KSP_CG_SYMMETRIC);
    KSPSetInitialGuessNonzero(solver,PETSC_TRUE);
    KSPSetUp(solver);
  }
  int solve(const double *b, double *x, size_t rhs, boost::property_tree::ptree &opts) {
    KSPSetTolerances(solver,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    vec_create_seq_with_array(comm, 1, Dim, b, &B_);
    vec_create_seq_with_array(comm, 1, Dim, x, &X_);
    KSPSolve(solver, B_, X_);

    vec_destroy(&B_);
    vec_destroy(&X_);
  }
  ~PETsc_CG_imp() {
    // std::cout<< "PETsc_CG_imp finalized" << std::endl;
    mat_destroy(&A_);
    KSP_destroy(&solver);
  }
private:
  Mat A_;
  Vec X_, B_;
  zjucad::matrix::matrix<int> ptr_, idx_;
  KSP solver;
  MPI_Comm comm;
  int Dim;
};

extern "C" {
PETsc* create_PETsc()
{
  return new PETsc_imp();
}
PETsc_CG* create_PETsc_CG(const double * val, const int32_t * idx,
			  const int32_t * ptr, const size_t nnz,
			  const size_t row, const size_t col,
			  const char *pchar)
{
  return new PETsc_CG_imp(val,idx, ptr, nnz,row, col, pchar);
}
}
