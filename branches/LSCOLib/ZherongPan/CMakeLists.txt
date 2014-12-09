cmake_minimum_required(VERSION 2.8)

project(LSCOLib)

file(GLOB header *.h)
file(GLOB source *.c *.cpp)

set(3rd_include "3rd part libraries include directories" CACHE PATH "3rd part libraries include directories") 
set(3rd_lib_path "3rd part libraries lib path" CACHE PATH "3rd part libraries lib path") 
set(eigen_path "eigen path" CACHE PATH "eigen Path")
set(CommonFilePZR_include "CommonFilePZR_include" CACHE PATH "CommonFilePZR_include")
set(CMAKE_CXX_FLAGS "-std=c++0x -fpermissive ${CMAKE_CXX_FLAGS}")

if(MSVC)
	if(NOT DEFINED 3rd_include)
		message(SEND_ERROR "need 3rd_include")
	endif(NOT DEFINED 3rd_include)

	if(NOT DEFINED eigen_path)
		message(SEND_ERROR "need eigen_path")
	endif(NOT DEFINED eigen_path)

	if(NOT DEFINED CommonFilePZR_include)
		message(SEND_ERROR "need CommonFilePZR")
	endif(NOT DEFINED CommonFilePZR_include)

	include_directories(${3rd_include})
	include_directories(${eigen_path})
	include_directories(${CommonFilePZR_include})
else(MSVC)
	include_directories($ENV{HOME}/usr/include)
	include_directories($ENV{HOME}/usr/include/eigen3)
	include_directories($ENV{HOME}/usr/include/CommonFilePZR)
        include_directories(/usr/include/suitesparse)
endif(MSVC)
add_executable(test ${header} ${source} )

find_package(Boost 1.46 COMPONENTS system)
if(Boost_FOUND)
	target_link_libraries(test
                $ENV{HOME}/usr/lib/libipopt.so
                $ENV{HOME}/usr/lib/libcoinhsl.so
                $ENV{HOME}/usr/lib/libcoinmetis.so
                $ENV{HOME}/usr/lib/libcoinmumps.so
                gomp
                dl
                lapack
                blas
                $ENV{HOME}/usr/Linux/64/gcc4/lib/libhj-math.so
                $ENV{HOME}/usr/Linux/64/lib/libjtf-function.so
                $ENV{HOME}/usr/Linux/64/gcc4/lib/libjtf-mesh.so
                ${Boost_LIBRARIES}
 )
	message("${Boost_LIBRARIES}")
else(Boost_FOUND)
	message(SEND_ERROR "Need Boost.")
endif(Boost_FOUND)

install(DIRECTORY . 
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/LSCOLib
	FILES_MATCHING PATTERN "*.h" REGEX ".SVN.*" EXCLUDE)
