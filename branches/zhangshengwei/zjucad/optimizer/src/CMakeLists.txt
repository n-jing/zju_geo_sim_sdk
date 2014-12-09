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
        $ENV{HOME}/usr/include
  )
endif(WIN32)

add_library(zjucad-optimizer SHARED ${source})

target_link_libraries(zjucad-optimizer
  gmodule-2.0
  hj-sparse-solver
  hj-sparse-util
  zjucad-ptree
	zjucad-linear-solver
  lapack
  alglib
  libLBFGS
  HLBFGS-static
)

geo_sim_sdk_install_cxx(TARGETS zjucad-optimizer)

