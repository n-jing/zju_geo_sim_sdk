cmake_minimum_required(VERSION 2.8)

project(hj_geometry)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_CXX_FLAGS "-std=c++0x")

add_subdirectory(src)

geo_sim_sdk_install_header(hjlib/geometry)
