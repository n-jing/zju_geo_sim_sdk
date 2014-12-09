cmake_minimum_required(VERSION 2.8)
project(tri_face)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
INCLUDE_DIRECTORIES($ENV{HOME}/usr/include)
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
set(CMAKE_CXX_FLAGS "-fpermissive")
find_package( Boost 1.42.0 COMPONENTS date_time filesystem system program_options )

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message(${Boost_LIBRARIES})
else()
  message("--[info] can not find Boost")
endif()

FILE(GLOB HEADFILE *.h)
FILE(GLOB SOURCEFILE *.cpp)
add_subdirectory(test)
ADD_LIBRARY(zjucad-tri_sub-static STATIC ${HEADFILE} ${SOURCEFILE})
ADD_LIBRARY(zjucad-tri_sub SHARED ${HEADFILE} ${SOURCEFILE})
target_link_libraries(zjucad-tri_sub
  ${Boost_LIBRARIES})


target_link_libraries(zjucad-tri_sub
  ${Boost_LIBRARIES})

install(FILES "tri_subdivision.h"
  DESTINATION include/sxxlib/tri_sub)

geo_sim_sdk_install_cxx(TARGETS zjucad-tri_sub)
geo_sim_sdk_install_cxx(TARGETS zjucad-tri_sub-static)



