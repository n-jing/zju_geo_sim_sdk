add_executable(linear_solver_test main.cpp)
target_link_libraries( linear_solver_test
  hj-sparse-solver
  hj-sparse-util
  zjucad-ptree
  zjucad-linear-solver
  lapack
  alglib
  libLBFGS
  # petsc petscmat petscvec petscksp petscdm
  HLBFGS-static
  # mpi_cxx
)