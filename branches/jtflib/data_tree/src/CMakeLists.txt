aux_source_directory(. sources)

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

add_library(hj_data_tree STATIC ${sources})

geo_sim_sdk_install_cxx(TARGETS hj_data_tree)
