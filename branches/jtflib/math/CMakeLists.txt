cmake_minimum_required(VERSION 2.8)

project(jtf-math)

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
add_library(jtf-math SHARED ${header} ${source})
#add_executable(jtf-math  ${header} ${source})

# TODO: why this program can only compiled as linux 64 gcc4?
install(TARGETS jtf-math
  DESTINATION Linux/64/gcc4/lib/)

install(DIRECTORY include/
  DESTINATION include/jtflib/math
  PATTERN .svn EXCLUDE)

