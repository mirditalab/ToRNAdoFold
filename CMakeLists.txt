cmake_minimum_required(VERSION 3.14)
project(simdmfe LANGUAGES CXX)

# SimdMFE is a header-only Turner 2004 SIMD MFE folder (everything lives in
# src/*.h). Consumers embed it with:
#
#   add_subdirectory(path/to/rna EXCLUDE_FROM_ALL)   # headers only, no binaries
#   target_link_libraries(myapp PRIVATE simdmfe::simdmfe)
#
# EXCLUDE_FROM_ALL keeps the CLI drivers out of the parent's default build; the
# INTERFACE library has no build artifact of its own, so nothing compiles until
# a consumer actually #includes and uses the headers.

add_library(simdmfe INTERFACE)
add_library(simdmfe::simdmfe ALIAS simdmfe)
target_include_directories(simdmfe INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>)
target_compile_features(simdmfe INTERFACE cxx_std_11)

# Build the CLI drivers only when SimdMFE is the top-level project. Included via
# add_subdirectory() it stays header-only regardless of EXCLUDE_FROM_ALL; set
# -DSIMDMFE_BUILD_TOOLS=ON to build them as a subproject anyway.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(_simdmfe_tools_default ON)
else()
    set(_simdmfe_tools_default OFF)
endif()
option(SIMDMFE_BUILD_TOOLS "Build the SimdMFE CLI and verify binaries" ${_simdmfe_tools_default})

if(SIMDMFE_BUILD_TOOLS)
    add_executable(simdmfe_cli src/main.cpp)
    set_target_properties(simdmfe_cli PROPERTIES OUTPUT_NAME simdmfe)
    target_link_libraries(simdmfe_cli PRIVATE simdmfe::simdmfe)

    add_executable(simdmfe_verify src/verify.cpp)
    set_target_properties(simdmfe_verify PROPERTIES OUTPUT_NAME verify)
    target_link_libraries(simdmfe_verify PRIVATE simdmfe::simdmfe)

    # NEON-disabled build (scalar constant-factor baseline).
    add_executable(simdmfe_noneon src/main.cpp)
    target_compile_definitions(simdmfe_noneon PRIVATE DISABLE_NEON)
    target_link_libraries(simdmfe_noneon PRIVATE simdmfe::simdmfe)

    # Optional OpenMP batch build (parallel across sequences).
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        add_executable(simdmfe_omp src/main.cpp)
        target_link_libraries(simdmfe_omp PRIVATE simdmfe::simdmfe OpenMP::OpenMP_CXX)
    endif()
endif()
