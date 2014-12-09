aux_source_directory(. sources)

add_library(hj-math SHARED ${sources})

target_link_libraries(hj-math
  lapack
  qp
  blas
)

geo_sim_sdk_install_cxx(TARGETS hj-math)
