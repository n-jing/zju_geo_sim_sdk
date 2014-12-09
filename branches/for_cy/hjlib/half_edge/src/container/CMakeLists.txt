add_definitions(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

aux_source_directory(. source)

add_executable(half_edge_container ${source})

target_link_libraries(half_edge_container
	${Boost_LIBRARIES}
)
