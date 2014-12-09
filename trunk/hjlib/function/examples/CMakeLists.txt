if(MSVC)
  add_definitions("-DHJ_FUNCTION_API=__declspec(dllimport)")
endif(MSVC)

add_executable(hj-function-ex main.cpp)

target_link_libraries(hj-function-ex
  hj-function
  )
