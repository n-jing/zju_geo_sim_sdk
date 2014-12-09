cmake_minimum_required(VERSION 2.8)

project(ann-parallel)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(static)
add_subdirectory(shared)

install(DIRECTORY include/ANN
  DESTINATION ${OS}/include/
  PATTERN ".svn" EXCLUDE
  )
