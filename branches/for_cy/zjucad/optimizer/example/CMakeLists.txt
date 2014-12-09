file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(exp ${header} ${source})
target_link_libraries(exp 
	zjucad-optimizer)