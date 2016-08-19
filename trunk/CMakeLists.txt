cmake_minimum_required(VERSION 2.8)

include_directories(./)
include_directories(../3rd/official/ann/include/)

include_directories(${HOME}/include/laspack/)

link_directories(${CMAKE_SOURCE_DIR}/../3rd/official/ann/lib/)

add_subdirectory(hjlib)
add_subdirectory(jtflib)
add_subdirectory(zjucad)
