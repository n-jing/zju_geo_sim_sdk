aux_source_directory_with_headers(. sources headers)

add_library(csc_cholmod SHARED ${sources} ${headers})

target_link_libraries(csc_cholmod
  cholmod
)

set_target_properties(csc_cholmod PROPERTIES
  SOVERSION 2
  VERSION 2.205)

geo_sim_sdk_install_c(TARGETS csc_cholmod
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
