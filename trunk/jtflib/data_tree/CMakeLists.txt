cmake_minimum_required(VERSION 2.8)

project(data_tree)

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    include($ENV{HOME}\usr\share\cmake\Modules\geo_sim_sdk.cmake)
 elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if(MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
elseif(GCC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif(MSVC)

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
