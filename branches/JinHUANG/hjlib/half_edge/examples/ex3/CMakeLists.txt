aux_source_directory(. source)

include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_CXX_FLAGS "-std=c++0x")
add_executable(half_edge_ex3 ${source})

target_link_libraries(half_edge_ex3
	jtf-mesh )
