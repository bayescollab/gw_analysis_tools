cmake_minimum_required (VERSION 3.11)
project(gw_analysis_tools
	VERSION 0.0 	
	DESCRIPTION "Tools to help with gravitational wave analysis and prediction"
	LANGUAGES CXX C)

option(TEST_PROGRAMS "test-programs" OFF)
option(DOCS "DOCS" OFF)
option(IMPLEMENT_HDF5 "Implement HDF5 if available" ON)
option(DEBUG "set debugger options" OFF)
option(ENABLE_PV3_INIT_ERROR
	"Allow for IMRPhenomPv3 initialization to throw error. STDERR will report the error, which can slow down inference."
	OFF)

#SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
#SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++ -lc++abi")
if(DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall ")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -stdlib=libstdc++")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install/" CACHE PATH "..." FORCE)
endif()
# The version number.
set (GWAT_VERSION_MAJOR 1)
set (GWAT_VERSION_MINOR 0)

set (GWAT_ROOT_DIRECTORY "${PROJECT_SOURCE_DIR}/")
set (GWAT_INSTALL_DIRECTORY "${CMAKE_INSTALL_PREFIX}/")
message(STATUS "Install directory: ${GWAT_INSTALL_DIRECTORY}")

 
# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/gwat/GWATConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/gwat/GWATConfig.h"
  )

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	#Only compile if this project is the main project
	#Just leave empty for now
endif()
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/gwat/*.h")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

include_directories("${PROJECT_BINARY_DIR}/include/gwat")
#Required packages
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND AND IMPLEMENT_HDF5 )
	add_compile_definitions(_HDF5)
else()
	message ( "HDF5 is currently required for some MCMC and RJMCMC functionality" )
endif()

find_package(ADOLC REQUIRED)
find_package(OpenMP REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Armadillo REQUIRED)

###############################################################

find_package(BayesShip REQUIRED)
###############################################################

set(INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include/gwat;${PROJECT_SOURCE_DIR}/include;${OpenMP_INCLUDE_DIR};${ADOLC_INCLUDE_DIRS};${HDF5_INCLUDE_DIRS};${FFTW3_INCLUDE_DIR};${GSL_INCLUDE_DIR};${ARMADILLO_INCLUDE_DIRS}" )
set(LIB_DIRS "${HDF5_LIBRARIES};OpenMP::OpenMP_CXX;${ADOLC_LIBRARIES};${FFTW3_LIBRARIES};${GSL_LIBRARIES};m;${ARMADILLO_LIBRARIES}")

#Library
add_subdirectory("${CMAKE_SOURCE_DIR}/src")

#Executables
add_subdirectory("${CMAKE_SOURCE_DIR}/apps")

#Python
#add_subdirectory("${CMAKE_SOURCE_DIR}/gwatpy")

if(TEST_PROGRAMS)
	add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
endif()


find_package(Doxygen)
if(DOCS AND DOXYGEN_FOUND)
	#add_subdirectory(docs)	
	#add_subdirectory(docs_doxygen)	
endif()
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
#include_directories("${PROJECT_BINARY_DIR}")


add_subdirectory("${CMAKE_SOURCE_DIR}/data")

install(FILES ${HEADER_LIST} DESTINATION "include/gwat")
install(FILES "${CMAKE_BINARY_DIR}/include/gwat/GWATConfig.h" DESTINATION "include/gwat")
