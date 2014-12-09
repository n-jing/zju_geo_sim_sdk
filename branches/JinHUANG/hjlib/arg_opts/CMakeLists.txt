cmake_minimum_required(VERSION 2.8)

project(arg_opts)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

add_subdirectory(src)

geo_sim_sdk_install_header(hjlib/arg_opts)
