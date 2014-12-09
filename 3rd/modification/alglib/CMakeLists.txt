cmake_minimum_required(VERSION 2.8)

project(alglib)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

add_subdirectory(src)
add_subdirectory(tests)
