file(GLOB header *.h)
file(GLOB source *.cpp)
find_package(MPI)
add_executable(example ${header} ${source})
target_link_libraries(example jtf-optimizer
zjucad-ptree
zjucad-optimizer
zjucad-linear-solver
${MPI_CXX_LIBRARIES}
)
