cmake_minimum_required(VERSION 2.8)

project(zjucad-ptree)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
  message("--[info] can not find Boost")
endif()


file(GLOB HEADER include/*.h)
file(GLOB SOURCE src/*.cpp)

include_directories(include)

add_library(zjucad-ptree SHARED ${SOURCE} ${HEADER})
add_library(zjucad-ptree-static STATIC ${SOURCE} ${HEADER})

target_link_libraries(zjucad-ptree)

install(DIRECTORY include/
  DESTINATION include/zjucad/ptree
  PATTERN .svn EXCLUDE)
geo_sim_sdk_install_cxx(TARGETS zjucad-ptree)
geo_sim_sdk_install_cxx(TARGETS zjucad-ptree-static)
