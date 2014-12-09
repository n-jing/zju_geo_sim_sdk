aux_source_directory(. source)
file(GLOB header *.h)

if(WIN32)
else(WIN32)
  include_directories(
	/usr/include/petsc
	/usr/lib/openmpi/include
	/usr/lib/openmpi/lib
        /usr/include/glib-2.0
        $ENV{HOME}/usr/include
  )
endif(WIN32)

add_library(zjucad-optimizer SHARED ${source} ${header})

target_link_libraries(zjucad-optimizer
  hj-sparse-util
  zjucad-ptree
	zjucad-linear-solver
  lapack
  alglib
  libLBFGS
  HLBFGS-static
)

geo_sim_sdk_install_cxx(TARGETS zjucad-optimizer)

