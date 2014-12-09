aux_source_directory(. sources)

if(WIN32)
  add_definitions("-DHJ_FUNCTION_API=__declspec(dllexport)")
endif(WIN32)

add_library(hj-function SHARED ${sources})

geo_sim_sdk_install_c(TARGETS hj-function)
