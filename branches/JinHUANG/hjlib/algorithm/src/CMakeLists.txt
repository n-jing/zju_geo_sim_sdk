aux_source_directory(. sources)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(hj_heap ${sources})
