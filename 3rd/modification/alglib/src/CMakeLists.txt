include_directories(${CMAKE_CURRENT_SOURCE_DIR})

aux_source_directory(. sources)

add_library(alglib SHARED ${sources})

set_target_properties(alglib PROPERTIES
  SOVERSION 3
  VERSION 3.4.0)

file(GLOB hdrs ./ *.h)

geo_sim_sdk_install_cxx(TARGETS alglib)

install(FILES ${hdrs}
  DESTINATION include)
