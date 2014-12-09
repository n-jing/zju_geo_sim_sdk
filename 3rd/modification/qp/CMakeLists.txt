cmake_minimum_required(VERSION 2.8)

project(qp)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

set(qp-source
  qld.c
  qp.c
)

add_library(qp-static STATIC ${qp-source})
add_library(qp SHARED ${qp-source})
target_link_libraries(qp
  f2c
)

add_executable(qp-example example.c)
target_link_libraries(qp-example
  qp)

geo_sim_sdk_install_c(TARGETS qp qp-static
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

install(FILES qp.h
  DESTINATION include/
)





