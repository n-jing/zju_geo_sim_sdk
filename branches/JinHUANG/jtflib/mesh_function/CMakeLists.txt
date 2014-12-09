cmake_minimum_required(VERSION 2.8)

project(mesh_function)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

file(GLOB header *.h)
add_custom_target(mesh_function SOURCES ${header})
install(FILES ${sources}
  DESTINATION include/jtflib/mesh_function
  )
