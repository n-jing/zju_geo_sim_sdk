aux_source_directory_with_headers(. sources headers)

add_library(csc_umfpack SHARED ${sources} ${headers})

target_link_libraries(csc_umfpack
   ${UMFPACK_LIBRARIES}
)

set_target_properties(csc_umfpack PROPERTIES
  SOVERSION 2
  VERSION 2.205)

geo_sim_sdk_install_c(TARGETS csc_umfpack
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
