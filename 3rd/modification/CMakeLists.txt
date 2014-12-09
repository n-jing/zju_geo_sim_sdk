cmake_minimum_required(VERSION 2.8)

project(3rd-lib)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

add_subdirectory(alglib EXCLUDE_FROM_ALL)
add_subdirectory(ann EXCLUDE_FROM_ALL)
add_subdirectory(HLBFGS EXCLUDE_FROM_ALL)
add_subdirectory(laspack EXCLUDE_FROM_ALL)
add_subdirectory(liblbfgs EXCLUDE_FROM_ALL)
add_subdirectory(qp EXCLUDE_FROM_ALL)
add_subdirectory(mesquite-2.2.0 EXCLUDE_FROM_ALL)
add_subdirectory(verdict EXCLUDE_FROM_ALL)
