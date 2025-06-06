
file(GLOB SRC_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")

#set(INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include/gwat;${OpenMP_INCLUDE_DIR};${ADOLC_INCLUDE_DIR};${HDF5_INCLUDE_DIR};${FFTW3_INCLUDE_DIR};${GSL_INCLUDE_DIR};${PTRJMCMC_INCLUDE_DIR}" )
#set(LIB_DIRS "${HDF5_LIBRARIES};OpenMP::OpenMP_CXX;${ADOLC_LIBRARIES};${FFTW3_LIBRARIES};${GSL_LIBRARIES};m;${PTRJMCMC_LIBRARIES}")

add_library(OBJS OBJECT ${SRC_LIST})
set_property(TARGET OBJS PROPERTY  POSITION_INDEPENDENT_CODE ON)
set_property(TARGET OBJS PROPERTY CXX_STANDARD 11)
target_link_libraries( OBJS PUBLIC "${LIB_DIRS}")
target_include_directories( OBJS PUBLIC ${INCLUDE_DIRS})

add_library(gwat SHARED $<TARGET_OBJECTS:OBJS> ${HEADER_LIST})
add_library(gwat_static STATIC $<TARGET_OBJECTS:OBJS> ${HEADER_LIST})
set_property(TARGET gwat PROPERTY CXX_STANDARD 11)
set_property(TARGET gwat_static PROPERTY CXX_STANDARD 11)

target_link_libraries( gwat PUBLIC "${LIB_DIRS}")
target_include_directories( gwat PUBLIC ${INCLUDE_DIRS})
target_link_libraries( gwat_static PUBLIC "${LIB_DIRS}")
target_include_directories( gwat_static PUBLIC ${INCLUDE_DIRS})



#Build and name
set_target_properties(gwat_static PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set_target_properties(gwat PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set_target_properties(gwat_static PROPERTIES OUTPUT_NAME "gwat")

#Install
install(TARGETS gwat DESTINATION lib)
install(TARGETS gwat_static DESTINATION lib)
