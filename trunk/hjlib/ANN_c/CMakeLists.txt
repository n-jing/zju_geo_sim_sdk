cmake_minimum_required(VERSION 2.8)

project(hj_ANN_c)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()
link_geo_sim_sdk()

aux_source_directory(. sources)
add_library(hj_ANN_c SHARED ${sources})

target_link_libraries(hj_ANN_c ANN)

geo_sim_sdk_install_c(TARGETS hj_ANN_c)

install(FILES ANN_c.h
  DESTINATION include/hjlib)
