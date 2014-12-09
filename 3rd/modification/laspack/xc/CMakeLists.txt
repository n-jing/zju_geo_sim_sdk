aux_source_directory(. sources)

add_library(laspack-xc SHARED ${sources})
add_library(laspack-xc-static STATIC ${sources})
