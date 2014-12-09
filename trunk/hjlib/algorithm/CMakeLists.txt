cmake_minimum_required(VERSION 2.8)

project(hj_heap)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

add_subdirectory(src)

geo_sim_sdk_install_header(hjlib/algorithm)
