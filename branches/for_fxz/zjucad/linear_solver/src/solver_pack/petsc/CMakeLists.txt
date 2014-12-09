include_directories(/usr/include/petsc)

add_library(dl_petsc SHARED ../../dl_petsc.h dl_petsc.cpp)

target_link_libraries(dl_petsc
  petsc petscmat petscvec petscksp petscdm
)

geo_sim_sdk_install_cxx(TARGETS dl_petsc)
