include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_executable(test_c test_c.cpp)
target_link_libraries(test_c
  alglib)

# add_executable(test_i test_i.cpp)
# target_link_libraries(test_i
#  alglib)
