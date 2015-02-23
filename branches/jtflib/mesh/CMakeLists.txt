cmake_minimum_required(VERSION 2.8)

project(jtf-mesh)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    include($ENV{HOME}\usr\share\cmake\Modules\geo_sim_sdk.cmake)
 elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")


set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
include_geo_sim_sdk()
link_geo_sim_sdk()

file(GLOB header include/*.h)
file(GLOB source src/*.cpp)

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
    message("CAN NOT FIND BOOST")
endif(Boost_FOUND)
add_library(jtf-mesh SHARED ${header} ${source})
target_link_libraries(jtf-mesh jtf-util lapack)

install(TARGETS jtf-mesh
  DESTINATION ${CMAKE_INSTALL_PREFIX}/Linux/64/gcc4/lib/)

install(DIRECTORY include/
  DESTINATION include/jtflib/mesh
  PATTERN .svn EXCLUDE)

