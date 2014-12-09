project(mk_ptree)

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package( Boost 1.42.0 COMPONENTS date_time filesystem system program_options )

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message(${Boost_LIBRARIES})
else()
  message("--[info] can not find Boost")
endif()


file(GLOB HEADER ./*.h)
file(GLOB SOURCE ./*.cpp)
#add_library(jtf-ptree static ${SOURCE} ${HEADER})
add_executable(mk_ptree ${SOURCE} ${HEADER})
target_link_libraries(mk_ptree
  ${Boost_LIBRARIES})


install(TARGETS mk_ptree 
  DESTINATION ${CMAKE_INSTALL_PREFIX}/Linux/64/gcc4/bin/)
