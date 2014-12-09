cmake_minimum_required(VERSION 2.8)

project(tnt)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()

if(MSVC)
else(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS
    "-static-libgcc -s")
  set(CMAKE_SHARED_LINKER_FLAGS
    "-static-libgcc -s")
endif(MSVC)


install(DIRECTORY include/
  DESTINATION include/tnt
  FILES_MATCHING PATTERN "*.h" REGEX ".svn.*" EXCLUDE)
