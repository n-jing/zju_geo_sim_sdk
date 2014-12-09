cmake_minimum_required(VERSION 2.8)

project(zjucad-optimizer)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

add_subdirectory(src)

install(DIRECTORY include/
  DESTINATION include/zjucad/optimizer
  PATTERN .svn EXCLUDE)

