cmake_minimum_required(VERSION 2.6)

project(hj-sparse)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

if(MSVC)
else(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS
    "-static-libgcc -s")
  set(CMAKE_SHARED_LINKER_FLAGS
    "-static-libgcc -s")
endif(MSVC)

add_subdirectory(src)
add_subdirectory(examples)

geo_sim_sdk_install_header(hjlib/sparse)
