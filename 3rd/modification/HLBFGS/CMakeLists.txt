cmake_minimum_required(VERSION 2.8)

project(HLBFGS)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

if(MSVC)
else(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS
    "-static-libgcc -s")
  set(CMAKE_SHARED_LINKER_FLAGS
    "-static-libgcc -s")
endif(MSVC)

aux_source_directory(. sources)

add_library(HLBFGS-static STATIC ${sources})
add_library(HLBFGS SHARED ${sources})

file(GLOB headers "*.h")
install(FILES ${headers} DESTINATION "include/HLBFGS")
geo_sim_sdk_install_cxx(TARGETS HLBFGS-static HLBFGS)
