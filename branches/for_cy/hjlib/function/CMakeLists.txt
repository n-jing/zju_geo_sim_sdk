cmake_minimum_required(VERSION 2.8)

project(hj-function)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

include_geo_sim_sdk()

if(MSVC)
else(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS
    "-static-libgcc -s")
  set(CMAKE_SHARED_LINKER_FLAGS
    "-static-libgcc -s")
endif(MSVC)

add_subdirectory(src)
add_subdirectory(examples)

geo_sim_sdk_install_header(hjlib/function)
