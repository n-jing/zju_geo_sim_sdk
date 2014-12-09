file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(meshmodel SHARED ${HEADERS} ${SOURCES})

install(DIRECTORY ./
  DESTINATION ${CMAKE_SOURCE_DIR}/include/MeshModel
  FILES_MATCHING PATTERN "*.h" REGEX ".svn.*" EXCLUDE)
