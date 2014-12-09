file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(common SHARED ${HEADERS} ${SOURCES})

install(DIRECTORY ./
  DESTINATION ${CMAKE_SOURCE_DIR}/include/Common
  FILES_MATCHING PATTERN "*.h" REGEX ".svn.*" EXCLUDE)
