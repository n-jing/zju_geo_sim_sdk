cmake_minimum_required(VERSION 2.8)

project(libLBFGS)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

aux_source_directory(lib sources)

add_library(libLBFGS-static STATIC ${sources})
add_library(libLBFGS SHARED ${sources})

file(GLOB headers "include/*.h")
install(FILES ${headers} DESTINATION "include/libLBFGS")
geo_sim_sdk_install_c(TARGETS libLBFGS-static libLBFGS)
