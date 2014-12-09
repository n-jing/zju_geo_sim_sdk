cmake_minimum_required(VERSION 2.8)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
set(subdir util mesh mesh_function algorithm function optimizer math)
foreach(dir ${subdir})
add_subdirectory(${dir})
endforeach()

