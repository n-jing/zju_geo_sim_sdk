aux_source_directory(. source)

include_geo_sim_sdk()
link_geo_sim_sdk()

add_executable(half_edge_ex2 ${source})

target_link_libraries(half_edge_ex2
	jtf-mesh)
