aux_source_directory(. sources)

add_library(hj_data_tree SHARED ${sources})

geo_sim_sdk_install_cxx(TARGETS hj_data_tree)
