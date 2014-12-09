cmake_minimum_required(VERSION 2.8)

project(metis)

set(CMAKE_C_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0  -g -ggdb")
#jtf added -D__VC__
set(CMAKE_C_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -DNDEBUG -D__VC__") 

aux_source_directory(${CMAKE_SOURCE_DIR}/Lib metis-lib-source)
include_directories(${CMAKE_SOURCE_DIR}/Lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_library(metis SHARED ${metis-lib-source})
set_target_properties(metis PROPERTIES
  SOVERSION 4
  VERSION 4.0
)

add_library(metis-static STATIC ${metis-lib-source})
