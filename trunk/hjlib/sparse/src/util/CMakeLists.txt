aux_source_directory_with_headers(. sources headers)

add_library(hj-sparse-util SHARED ${sources} ${headers})

target_link_libraries(hj-sparse-util
  cholmod umfpack)

set_target_properties(hj-sparse-util PROPERTIES
  SOVERSION 2
  VERSION 2.205)

geo_sim_sdk_install_c(TARGETS hj-sparse-util
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
