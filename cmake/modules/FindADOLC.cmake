# - Try to find Adolc
# Once done this will define
#  ADOLC_FOUND - System has Adolc
#  ADOLC_INCLUDE_DIRS - The Adolc include directories
#  ADOLC_LIBRARY_DIRS - The library directories needed to use Adolc
#  ADOLC_LIBRARIES    - The libraries needed to use Adolc
#  Adapted from https://raw.githubusercontent.com/joaoleal/CppADCodeGen/master/cmake/FindADOLC.cmake

IF (ADOLC_INCLUDES AND ADOLC_LIBRARIES)
  SET(ADOLC_FIND_QUIETLY TRUE)
ENDIF (ADOLC_INCLUDES AND ADOLC_LIBRARIES)


FIND_PACKAGE(PkgConfig)
IF( PKG_CONFIG_FOUND )
  pkg_check_modules( ADOLC QUIET adolc>=2.3)
ENDIF( PKG_CONFIG_FOUND )


IF( ADOLC_FOUND )
  IF(NOT ADOLC_FIND_QUIETLY)
    MESSAGE(STATUS "package adolc ${ADOLC_VERSION} found")
  ENDIF()
ELSE( ADOLC_FOUND )
  string(REPLACE ":" ";" INC_SEARCH_LIST "$ENV{CPATH}")
  FIND_PATH(ADOLC_INCLUDE_DIRS NAMES adolc
            HINTS  
		   "$ENV{ADOLC_HOME}/include"
                   "/usr/include" 
		   ${INC_SEARCH_LIST}
		   )
           
  string(REPLACE ":" ";" LIB_SEARCH_LIST "$ENV{LIBRARY_PATH}")
  FIND_LIBRARY(ADOLC_LIBRARY adolc
                HINTS "$ENV{ADOLC_HOME}/lib64"
                      "/usr/lib" 
		      ${LIB_SEARCH_LIST})
  #CHECK_SYMBOL_EXISTS(function ADOLC_LIBRARY FUNCTION_FOUND)
   
  SET(ADOLC_INCLUDE_DIR ${ADOLC_INCLUDE_DIRS})
  SET(ADOLC_LIBRARIES ${ADOLC_LIBRARY})

  INCLUDE(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set LIBIPOPT_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(ADOLC  DEFAULT_MSG
                                    ADOLC_INCLUDE_DIRS ADOLC_LIBRARIES)

  MARK_AS_ADVANCED(ADOLC_INCLUDE_DIRS ADOLC_LIBRARIES)
  
  IF( ADOLC_FOUND AND NOT ADOLC_FIND_QUIETLY )
    MESSAGE(STATUS "package adolc found")
  ENDIF()
ENDIF( ADOLC_FOUND )
