aux_source_directory(. sources)

add_library(hj_math_func SHARED ${sources})

geo_sim_sdk_install_cxx(TARGETS hj_math_func)
