cmake_minimum_required(VERSION 2.8)

project(hjlib_half_edge)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()
link_geo_sim_sdk()
include_directories($ENV{HOME}/usr/include)

set(CMAKE_CXX_FLAGS "-fpermissive -std=c++0x")

add_subdirectory(doc) # documentation

add_subdirectory(src) # unit test

add_subdirectory(examples)



geo_sim_sdk_install_header(hjlib/half_edge)

