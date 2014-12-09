aux_source_directory(. sources)

add_executable(hj_data_tree_sample ${sources})
target_link_libraries(hj_data_tree_sample
hj_data_tree)
