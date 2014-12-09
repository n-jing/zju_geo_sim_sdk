add_executable(hj-sparse-vec vec.cpp)
add_executable(hj-sparse-example main.cpp)

if(WIN32)
else(WIN32)
target_link_libraries(hj-sparse-example
  gfortran
)
endif(WIN32)
