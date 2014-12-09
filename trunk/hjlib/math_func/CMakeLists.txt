cmake_minimum_required(VERSION 2.8)

project(hj_math_func)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

set(CMAKE_CXX_FLAGS "-std=c++0x -fPIC -fopenmp ${CMAKE_CXX_FLAGS}")

include_geo_sim_sdk()

add_subdirectory(src)

geo_sim_sdk_install_header(hjlib/math_func)
