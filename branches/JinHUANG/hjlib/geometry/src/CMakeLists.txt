aux_source_directory(. sources)

add_library(hj_geometry SHARED ${sources})

target_link_libraries(hj_geometry
	lapack)

geo_sim_sdk_install_cxx(TARGETS hj_geometry)
