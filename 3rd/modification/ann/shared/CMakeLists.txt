aux_source_directory(../src sources)

add_definitions(-DDLL_EXPORTS)
add_library(ann-parallel SHARED ${sources})

geo_sim_sdk_install_cxx(TARGETS ann-parallel)
