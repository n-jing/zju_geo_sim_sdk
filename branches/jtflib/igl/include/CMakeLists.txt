file(GLOB header *.h)
include_directories(/usr/include/eigen3)
add_custom_target(igl_include SOURCES ${header})
install(DIRECTORY .
	DESTINATION include/jtflib/igl
	PATTERN .svn EXCLUDE)
