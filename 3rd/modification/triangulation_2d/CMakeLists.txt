cmake_minimum_required(VERSION 2.8)

project(triangulation_2d)

set(CMAKE_CXX_FLAGS "-fPIC -fpermissive")

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

if(MSVC)
else(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS
    "-static-libgcc -s")
  set(CMAKE_SHARED_LINKER_FLAGS
    "-static-libgcc -s")
endif(MSVC)

aux_source_directory(. sources)

add_library(triangulation_2d-static STATIC ${sources})
add_library(triangulation_2d SHARED ${sources})

file(GLOB headers "*.h")
install(FILES ${headers} DESTINATION "include/triangulation_2d")
geo_sim_sdk_install_cxx(TARGETS triangulation_2d-static triangulation_2d)
