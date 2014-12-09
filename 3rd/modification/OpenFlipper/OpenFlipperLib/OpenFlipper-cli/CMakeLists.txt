include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/OpenMesh/src
  ${GLEW_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/OpenFlipperLib
)

aux_source_directory(${CMAKE_SOURCE_DIR}/OpenFlipper-cli sources)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/Build/bin)
add_executable(OpenFlipper-cli ${sources})

target_link_libraries(OpenFlipper-cli
  OpenFlipperLib
)
