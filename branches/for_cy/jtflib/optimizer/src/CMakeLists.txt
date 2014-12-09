file(GLOB header *.h)
file(GLOB source *.cpp)

add_library(jtf-optimizer SHARED ${source} ${header})

target_link_libraries( jtf-optimizer
  petsc petscmat petscvec petscksp petscdm
  mpi_cxx
  alglib
  jtf-function
	zjucad-ptree
	zjucad-optimizer
	zjucad-linear-solver
)

install(TARGETS jtf-optimizer
  DESTINATION ${CMAKE_INSTALL_PREFIX}/Linux/64/gcc4/lib/)
