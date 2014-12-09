cmake_minimum_required(VERSION 2.8)

project(zjucad-linear-solver)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
add_subdirectory(src)
#add_subdirectory(example)

install(DIRECTORY include/
  DESTINATION include/zjucad/linear_solver
  PATTERN .svn EXCLUDE)



