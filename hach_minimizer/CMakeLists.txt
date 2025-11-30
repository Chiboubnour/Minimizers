cmake_minimum_required(VERSION 3.15)
project(Minimizers LANGUAGES CXX)

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options globales Linux
add_compile_options(-Wall -Wextra -Werror -O3)

# March native si supporté
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
if (HAS_MARCH_NATIVE)
    add_compile_options(-march=native)
endif()

# Filesystem obligatoire
link_libraries(stdc++fs)

# OpenMP
find_package(OpenMP REQUIRED)

# ---- LIBRAIRIE INTERFACE commune ----
add_library(common_interface INTERFACE)
target_include_directories(common_interface INTERFACE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/hls
)

# ---- EXEC MONO ----
add_executable(minimizer_mono
    test/main_mono.cpp
    src/hach_minimizer_mono.cpp
)
target_link_libraries(minimizer_mono PRIVATE common_interface)

# ---- EXEC MULTI (OpenMP) ----
add_executable(minimizer_multi
    test/main_multi.cpp
    src/hach_minimizer_multi.cpp
)
target_link_libraries(minimizer_multi PRIVATE 
    common_interface
    OpenMP::OpenMP_CXX
)

# ---- EXEC HLS ----
add_executable(minimizer_hls
    test/main_hls.cpp
    hls/hach_minimizer_hls.cpp
)
target_link_libraries(minimizer_hls PRIVATE common_interface)
target_compile_options(minimizer_hls PRIVATE -Wno-sign-compare -Wno-conversion)

# ---- Génération des données ----
add_executable(generate_minimizers dataset/main_data.cpp)
target_link_libraries(generate_minimizers PRIVATE 
    common_interface
    OpenMP::OpenMP_CXX
)

# ---- TESTS ----
enable_testing()

add_custom_target(generate_test_data
    COMMAND generate_minimizers
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_dependencies(minimizer_mono generate_test_data)

add_test(NAME test_minimizer_mono
    COMMAND minimizer_mono ${CMAKE_BINARY_DIR}/generated/minimizers_32MB.bin 1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
