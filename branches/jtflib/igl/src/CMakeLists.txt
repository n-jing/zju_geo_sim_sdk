file(GLOB source *.cpp)
include_directories($ENV{HOME}/usr/include)
include_directories(/usr/include/eigen3)
include_directories($ENV{HOME}/usr/include/AntTweakBar)
add_library(jtf-igl SHARED ${source})
install(TARGETS jtf-igl
  DESTINATION ${CMAKE_INSTALL_PREFIX}/Linux/64/gcc4/lib/)
