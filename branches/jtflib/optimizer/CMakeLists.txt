cmake_minimum_required(VERSION 2.8)

project(jtf-optimizer)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    include($ENV{HOME}\usr\share\cmake\Modules\geo_sim_sdk.cmake)
 elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# TODO: I think the following only works on Linux
include_directories(
/usr/include/petsc/
/usr/lib/openmpi/include
/usr/lib/openmpi/lib)


add_subdirectory(src)
add_subdirectory(example)

install(DIRECTORY include/
  DESTINATION include/jtflib/optimizer
  PATTERN .svn EXCLUDE)

