cmake_minimum_required(VERSION 2.8)
project(jtf-util)

file(GLOB header include/*.h)
file(GLOB source src/*.cpp)


add_library(jtf-util SHARED ${header} ${source})

install(TARGETS jtf-util
 DESTINATION ${CMAKE_INSTALL_PREFIX}/Linux/64/gcc4/lib/) 
install(FILES ${header}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/jtflib/util )
