cmake_minimum_required(VERSION 2.8)

project(zjucad_matrix)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)


add_subdirectory(src)

file(GLOB headers include/*.h)

install(FILES ${headers}
  DESTINATION include/zjucad/matrix)
