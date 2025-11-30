cmake_minimum_required(VERSION 3.10)

project(Minimizers VERSION 1.0 LANGUAGES CXX)

# Options de configuration
option(ENABLE_HLS "Enable HLS-specific code" OFF)
option(BUILD_TESTS "Build test executables" ON)

# Configuration pour Linux
if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
    add_definitions(-DLINUX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

# Configuration du compilateur
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
    -Wall
    -Wextra
    -O2
    -std=c++11
    -I${CMAKE_SOURCE_DIR}/include
    -Wno-sign-compare
    -Wno-unused-parameter
    #-Wno-deprecated-copy-with-user-provided-copy
)


# Inclure les répertoires nécessaires
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/hls
    ${CMAKE_SOURCE_DIR}/src
)

# Définition des fichiers sources
set(COMMON_SOURCES
    ${CMAKE_SOURCE_DIR}/src/ref_minimizer.cpp
    ${CMAKE_SOURCE_DIR}/src/thr_minimizer_v1.cpp
    ${CMAKE_SOURCE_DIR}/hls/minimizer_v1.cpp
)

# Exécutable de référence
add_executable(minimizer_ref
    ${CMAKE_SOURCE_DIR}/test/main_ref.cpp
    ${CMAKE_SOURCE_DIR}/src/ref_minimizer.cpp
)

# Exécutable version threadée v1
add_executable(minimizer_thr_v1
    ${CMAKE_SOURCE_DIR}/test/main_thr_v1.cpp
    ${CMAKE_SOURCE_DIR}/src/thr_minimizer_v1.cpp
)

# Exécutable version threadée v2 (si disponible)
if(EXISTS "${CMAKE_SOURCE_DIR}/test/main_thr_v2.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/src/thr_minimizer_v2.cpp")
    add_executable(minimizer_thr_v2
        ${CMAKE_SOURCE_DIR}/test/main_thr_v2.cpp
        ${CMAKE_SOURCE_DIR}/src/thr_minimizer_v2.cpp
    )
endif()

# Exécutable HLS v1 (si disponible)
if(EXISTS "${CMAKE_SOURCE_DIR}/test/main_hls_v1.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/hls/minimizer_v1.cpp")
    add_executable(minimizer_hls_v1
        ${CMAKE_SOURCE_DIR}/test/main_hls_v1.cpp
        ${CMAKE_SOURCE_DIR}/hls/minimizer_v1.cpp
    )
    target_compile_definitions(minimizer_hls_v1 PRIVATE _PRAGMA_HLS_)
endif()

# Exécutable HLS v2 (si disponible)
if(EXISTS "${CMAKE_SOURCE_DIR}/test/main_hls_v2.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/hls/minimizer_v2.cpp")
    add_executable(minimizer_hls_v2
        ${CMAKE_SOURCE_DIR}/test/main_hls_v2.cpp
        ${CMAKE_SOURCE_DIR}/hls/minimizer_v2.cpp
    )
    target_compile_definitions(minimizer_hls_v2 PRIVATE _PRAGMA_HLS_)
endif()

# Options de liaison
foreach(target minimizer_ref minimizer_thr_v1 minimizer_thr_v2 minimizer_hls_v1 minimizer_hls_v2)
    if(TARGET ${target})
        target_link_libraries(${target} PRIVATE ${CMAKE_THREAD_LIBS_INIT})
        set_target_properties(${target} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
            ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        )
    endif()
endforeach()

# Installation
install(TARGETS minimizer_ref minimizer_thr_v1 DESTINATION bin)

# Tests
if(BUILD_TESTS AND EXISTS "${CMAKE_SOURCE_DIR}/test")
    enable_testing()
    
    # Test de base pour la version de référence
    add_test(NAME test_minimizer_ref
             COMMAND ${CMAKE_BINARY_DIR}/bin/minimizer_ref --input ${CMAKE_SOURCE_DIR}/dataset/b64/test_64.nuc -o ${CMAKE_BINARY_DIR}/test_ref_output.txt
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    
    # Test pour la version threadée
    add_test(NAME test_minimizer_thr_v1
             COMMAND ${CMAKE_BINARY_DIR}/bin/minimizer_thr_v1 --input ${CMAKE_SOURCE_DIR}/dataset/b64/test_64.nuc -o ${CMAKE_BINARY_DIR}/test_thr_output.txt
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    
    # Comparaison des sorties
    add_test(NAME compare_outputs
             COMMAND ${CMAKE_COMMAND} -E compare_files
             ${CMAKE_BINARY_DIR}/test_ref_output.txt
             ${CMAKE_BINARY_DIR}/test_thr_output.txt)
endif()

# Message de fin
message(STATUS "Configuration terminée. Pour compiler, exécutez 'make'")
message(STATUS "Pour installer: make install")
message(STATUS "Pour exécuter les tests: make test")
message(STATUS "Pour nettoyer: make clean")
