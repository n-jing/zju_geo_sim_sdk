cmake_minimum_required(VERSION 2.8)

project(jtf-algorithm)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()
link_geo_sim_sdk()

file(GLOB header include/*.h)

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
    message("CAN NOT FIND BOOST")
endif(Boost_FOUND)

aux_source_directory(. source)
add_subdirectory(src)
install(DIRECTORY include/
  DESTINATION include/jtflib/algorithm
  PATTERN .svn EXCLUDE)
