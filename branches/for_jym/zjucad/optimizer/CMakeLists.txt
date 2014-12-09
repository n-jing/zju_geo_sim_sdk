cmake_minimum_required(VERSION 2.8)

project(zjucad-optimizer)

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

add_subdirectory(src)
add_subdirectory(example)

install(DIRECTORY include/
  DESTINATION include/zjucad/optimizer
  PATTERN .svn EXCLUDE)

