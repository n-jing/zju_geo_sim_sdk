file(GLOB header *.h)
file(GLOB source *.cpp)

add_library(dl_base SHARED ${header} ${source})

target_link_libraries(dl_base
  cholmod umfpack
  laspack
)

geo_sim_sdk_install_cxx(TARGETS dl_base)
