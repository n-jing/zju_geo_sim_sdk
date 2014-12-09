cmake_minimum_required(VERSION 2.8)

project(laspack)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)

add_subdirectory(laspack)
add_subdirectory(xc)
