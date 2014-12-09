cmake_minimum_required(VERSION 2.8)

project(data_tree)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

set(CMAKE_CXX_FLAGS "-std=c++0x")

include_geo_sim_sdk()
link_geo_sim_sdk()

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
    message("CAN NOT FIND BOOST")
endif(Boost_FOUND)

add_subdirectory(src)
add_subdirectory(sample)

geo_sim_sdk_install_header(hjlib/data_tree)
