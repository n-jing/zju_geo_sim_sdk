aux_source_directory(. sources)

add_library(laspack SHARED ${sources})
add_library(laspack-static STATIC ${sources})

if(UNIX)
  set_target_properties(laspack-static
	PROPERTIES COMPILE_FLAGS -fPIC)
endif(UNIX)

add_subdirectory(examples)

geo_sim_sdk_install_c(TARGETS laspack laspack-static
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

file(GLOB hdr ./ *.h)
install(FILES ${hdr}
  DESTINATION include/laspack)
