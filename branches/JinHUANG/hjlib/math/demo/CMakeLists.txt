aux_source_directory(. sources)

add_executable(hj-math-demo ${sources})

target_link_libraries(hj-math-demo
  hj-math)
