aux_source_directory(../src sources)

add_library(ann-parallel-static STATIC ${sources})

geo_sim_sdk_install_cxx(TARGETS ann-parallel-static)
