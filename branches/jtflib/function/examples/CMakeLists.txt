aux_source_directory(. sources)

if(MSVC)
  add_definitions("-DHJ_FUNCTION_API=__declspec(dllimport)")
endif(MSVC)

add_executable(jtf-function-ex ${sources})

target_link_libraries(jtf-function-ex
  jtf-function
  )
