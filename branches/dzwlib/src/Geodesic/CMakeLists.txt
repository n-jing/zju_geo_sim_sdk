file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(geodesic SHARED ${HEADERS} ${SOURCES})

install(DIRECTORY ./
  DESTINATION ${CMAKE_SOURCE_DIR}/include/Geodesic
  FILES_MATCHING PATTERN "*.h" REGEX ".svn.*" EXCLUDE)
