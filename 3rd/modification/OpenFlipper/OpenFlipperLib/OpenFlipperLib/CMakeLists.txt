include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/OpenMesh/src
  ${GLEW_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
)

aux_source_directory(${CMAKE_SOURCE_DIR}/OpenFlipperLib sources)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/Build/bin)
add_library(OpenFlipperLib SHARED ${sources})

target_link_libraries(OpenFlipperLib
  PluginLibStatic
  ACGStatic
  IsotropicRemesher
  ${GLUT_LIBRARIES}
  ${GLEW_LIBRARY}
  ${QT_LIBRARIES}
)

if(WIN32)
  target_link_libraries(OpenFlipperLib
	OpenMeshCore
	)
else()
  target_link_libraries(OpenFlipperLib
	OpenMeshCoreStatic
	)
endif()
