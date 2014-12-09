aux_source_directory(. sources)

add_library(hj-arg_opts STATIC ${sources})

geo_sim_sdk_install_cxx(TARGETS hj-arg_opts)
