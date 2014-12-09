include_directories(/usr/include/petsc)

aux_source_directory(. source)

if(WIN32)
# include_directories("C:/MinGW/MinGW32/msys/1.0/home/wegatron/usr/mpich3.1/include")
# link_directories("C:/MinGW/MinGW32/msys/1.0/home/wegatron/usr/mpich3.1/lib")
else(WIN32)
  include_directories(
	/usr/lib/openmpi/include
	/usr/lib/openmpi/lib
  /usr/include/glib-2.0
	/usr/lib/x86_64-linux-gnu/glib-2.0/include/
  $ENV{HOME}/usr/include
  )
endif(WIN32)

add_library(zjucad-linear-solver SHARED ${source})

target_link_libraries(zjucad-linear-solver
  dl_base
  gmodule-2.0
  hj-sparse-util
  zjucad-ptree
  lapack
  alglib
  libLBFGS
  HLBFGS-static
)

geo_sim_sdk_install_cxx(TARGETS zjucad-linear-solver)

if(NOT WIN32)
	add_subdirectory(solver_pack)
endif(NOT WIN32)
