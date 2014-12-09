if(WIN32)
  set(mathlib)
else(WIN32)
  set(mathlib m)
endif(WIN32)

aux_source_directory(lastest lastest-sources)
add_executable(lastest ${lastest-sources})
target_link_libraries(lastest laspack-static laspack-xc-static ${mathlib})

aux_source_directory(matropt matropt-sources)
add_executable(matropt ${matropt-sources})
target_link_libraries(matropt laspack-static laspack-xc-static ${mathlib})

aux_source_directory(mlstest mlstest-sources)
add_executable(mlstest ${mlstest-sources})
target_link_libraries(mlstest laspack-static laspack-xc-static ${mathlib})

aux_source_directory(vectopt vectopt-sources)
add_executable(vectopt ${vectopt-sources})
target_link_libraries(vectopt laspack-static laspack-xc-static ${mathlib})
