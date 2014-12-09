find_package(Doxygen)
if (NOT DOXYGEN_FOUND)
  message(FATAL_ERROR 
    "Doxygen is needed to build the documentation. Please install it correctly")
endif()

configure_file(Doxyfile
  ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)

add_custom_target(half_edge_doc ALL 
	COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
	SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
