#include <iostream>
#include <petscksp.h> 
#include <private/matimpl.h>
#include <petscversion.h>
#include "../../dl_petsc.h"
#include "petsc_ver.h"
#include <hjlib/sparse/sparse.h>

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
  PETsc_CG_imp(const hj::sparse::csc<double, int32_t> &A, const char *pc_str) {
    std::cout << "call PETsc_CG_imp" << std::endl;
    comm = MPI_COMM_SELF;
    Dim = A.size(1);
    if(sizeof(int32_t) == sizeof(int)) {
        MatCreateSeqAIJWithArrays(
              comm,
              A.size(1), A.size(2),
              const_cast<int *>(&A.ptr()[0]), const_cast<int *>(&A.idx()[0]), const_cast<double *>(&A.val()[0]), &A_);
      }
    else {
        ptr_ = A.ptr();
        idx_ = A.idx();
        MatCreateSeqAIJWithArrays(
              comm,
              A.size(1), A.size(2), &ptr_[0], &idx_[0], const_cast<double *>(&A.val()[0]), &A_);
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
    VecCreateSeqWithArray(comm, Dim, b, &B_);
    VecCreateSeqWithArray(comm, Dim, x, &X_);
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
PETsc_CG* create_PETsc_CG(hj::sparse::csc<double, int32_t>& A, const char *pchar)
{
  return new PETsc_CG_imp(A, pchar);
}
}
