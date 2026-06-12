add_library(efpe STATIC src/pyefpe.cpp)
target_include_directories(efpe PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(efpe PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)
target_compile_options(efpe PRIVATE -O3)

install(TARGETS efpe ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

if(TEST_PROGRAMS)
    add_executable(efpe_smoke tests/smoke.cpp)
    set_property(TARGET efpe_smoke PROPERTY
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/bin")
    target_link_libraries(efpe_smoke PRIVATE efpe)
    set_target_properties(efpe_smoke PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON)
endif()
