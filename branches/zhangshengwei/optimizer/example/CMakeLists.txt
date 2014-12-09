file(GLOB header *.h)
file(GLOB source *.cpp)
find_package(MPI)

add_executable(example ${header} ${source})
target_link_libraries(example jtf-optimizer
zjucad-ptree
zjucad-optimizer
  petsc petscmat petscvec petscksp petscdm
  zjucad-linear-solver
  dl_base
  ${MPI_CXX_LIBRARIES}
)
