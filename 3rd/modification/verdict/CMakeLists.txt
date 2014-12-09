cmake_minimum_required(VERSION 2.8)

project(visit_verdict)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

file(GLOB source *.C)
file(GLOB header *.h)

add_definitions(-DVERDICT_EXPORTS)

add_library(visit_verdict ${source} ${header})

geo_sim_sdk_install_cxx(TARGETS visit_verdict)

install(FILES ${header}
  DESTINATION include/verdict)
#VISIT_INSTALL_TARGETS(visit_verdict)
#SET_TARGET_PROPERTIES(verdict PROPERTIES FOLDER thirdparty_builtin/visit_verdict)
