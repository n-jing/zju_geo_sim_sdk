cmake_minimum_required(VERSION 2.8)

project(jtf-function)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    include($ENV{HOME}\usr\share\cmake\Modules\geo_sim_sdk.cmake)
 elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")


include_geo_sim_sdk()

if(MSVC)
else(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS
    "-static-libgcc -s")
  set(CMAKE_SHARED_LINKER_FLAGS
    "-static-libgcc -s")
endif(MSVC)

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

add_subdirectory(include)
add_subdirectory(src)
#add_subdirectory(examples)


install(DIRECTORY include/
  DESTINATION include/jtflib/function
  PATTERN .svn EXCLUDE)
