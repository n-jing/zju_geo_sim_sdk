cmake_minimum_required(VERSION 2.8)

project(hj_vtk)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

file(GLOB headers include/*.h)

install(FILES ${headers}
  DESTINATION include/hjlib/vtk)
