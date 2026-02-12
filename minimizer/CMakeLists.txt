cmake_minimum_required(VERSION 3.10)

project(Minimizers VERSION 1.0 LANGUAGES CXX)

# Options
option(ENABLE_HLS "Enable HLS-specific code" OFF)
option(BUILD_TESTS "Build test executables" ON)

# Configuration Linux
if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
    add_definitions(-DLINUX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
    -Wall
    -Wextra
    -O2
    -std=c++17
    -I${CMAKE_SOURCE_DIR}/include
    -Wno-sign-compare
    -Wno-unused-parameter
    -Wno-int-in-bool-context
    -Wno-maybe-uninitialized
    -Wno-uninitialized
)
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/hls
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/../tools/src/files  # pour stream_gz_reader
)

# Bibliothèque zlib
find_package(ZLIB REQUIRED)

# --------------------------------------------------------------------
# Fichier source pour stream_gz_reader
# --------------------------------------------------------------------
set(STREAM_GZ_READER_SRC
    ${CMAKE_SOURCE_DIR}/../tools/src/files/gz/reader/stream_gz_reader.cpp
)

# --------------------------------------------------------------------
# Executables
# --------------------------------------------------------------------

# minimizer_ref
add_executable(minimizer_ref
    ${CMAKE_SOURCE_DIR}/test/main_ref.cpp
    ${CMAKE_SOURCE_DIR}/src/ref_minimizer.cpp
    ${STREAM_GZ_READER_SRC}   
)
target_link_libraries(minimizer_ref PRIVATE ${ZLIB_LIBRARIES})
target_include_directories(minimizer_ref PRIVATE ${ZLIB_INCLUDE_DIRS})

# minimizer_thr_v1
add_executable(minimizer_thr_v1
    ${CMAKE_SOURCE_DIR}/test/main_thr_v1.cpp
    ${CMAKE_SOURCE_DIR}/src/thr_minimizer_v1.cpp
    ${STREAM_GZ_READER_SRC}  
)
target_link_libraries(minimizer_thr_v1 PRIVATE ${ZLIB_LIBRARIES})
target_include_directories(minimizer_thr_v1 PRIVATE ${ZLIB_INCLUDE_DIRS})

# minimizer_thr_v2
if(EXISTS "${CMAKE_SOURCE_DIR}/test/main_thr_v2.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/src/thr_minimizer_v2.cpp")
    add_executable(minimizer_thr_v2
        ${CMAKE_SOURCE_DIR}/test/main_thr_v2.cpp
        ${CMAKE_SOURCE_DIR}/src/thr_minimizer_v2.cpp
        ${STREAM_GZ_READER_SRC} 
    )
    target_link_libraries(minimizer_thr_v2 PRIVATE ${ZLIB_LIBRARIES})
    target_include_directories(minimizer_thr_v2 PRIVATE ${ZLIB_INCLUDE_DIRS})
endif()

# minimizer_thr_v3
if(EXISTS "${CMAKE_SOURCE_DIR}/test/main_thr_v3.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/src/thr_minimizer_v3.cpp")
    add_executable(minimizer_thr_v3
        ${CMAKE_SOURCE_DIR}/test/main_thr_v3.cpp
        ${CMAKE_SOURCE_DIR}/src/thr_minimizer_v3.cpp
        ${STREAM_GZ_READER_SRC} 
    )
    target_link_libraries(minimizer_thr_v3 PRIVATE ${ZLIB_LIBRARIES})
    target_include_directories(minimizer_thr_v3 PRIVATE ${ZLIB_INCLUDE_DIRS})
endif()

# minimizer_thr_v3
if(EXISTS "${CMAKE_SOURCE_DIR}/test/test.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/src/thr_minimizer_v3.cpp")
    add_executable(test_thr_v3
        ${CMAKE_SOURCE_DIR}/test/test.cpp
        ${CMAKE_SOURCE_DIR}/src/thr_minimizer_v3.cpp
        ${CMAKE_SOURCE_DIR}/test/catch2v3/catch_amalgamated.cpp
    )
    target_link_libraries(test_thr_v3 PRIVATE ${ZLIB_LIBRARIES})
    target_include_directories(test_thr_v3 PRIVATE ${ZLIB_INCLUDE_DIRS})
endif()

# HLS v1
if(EXISTS "${CMAKE_SOURCE_DIR}/test/main_hls_v1.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/hls/minimizer_v1.cpp")
    add_executable(minimizer_hls_v1
        ${CMAKE_SOURCE_DIR}/test/main_hls_v1.cpp
        ${CMAKE_SOURCE_DIR}/hls/minimizer_v1.cpp
    )
    target_compile_definitions(minimizer_hls_v1 PRIVATE _PRAGMA_HLS_)
endif()

# HLS v2
if(EXISTS "${CMAKE_SOURCE_DIR}/test/main_hls_v2.cpp" AND EXISTS "${CMAKE_SOURCE_DIR}/hls/minimizer_v2.cpp")
    add_executable(minimizer_hls_v2
        ${CMAKE_SOURCE_DIR}/test/main_hls_v2.cpp
        ${CMAKE_SOURCE_DIR}/hls/minimizer_v2.cpp
    )
    target_compile_definitions(minimizer_hls_v2 PRIVATE _PRAGMA_HLS_)
endif()

foreach(target minimizer_ref minimizer_thr_v1 minimizer_thr_v2 minimizer_thr_v3 minimizer_hls_v1 minimizer_hls_v2)
    if(TARGET ${target})
        target_link_libraries(${target} PRIVATE ${CMAKE_THREAD_LIBS_INIT})
        set_target_properties(${target} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
            ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        )
    endif()
endforeach()

# --------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------
install(TARGETS minimizer_ref minimizer_thr_v1 minimizer_thr_v2 minimizer_thr_v3 DESTINATION bin)

# --------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------
if(BUILD_TESTS AND EXISTS "${CMAKE_SOURCE_DIR}/test")
    enable_testing()
    
    add_test(NAME test_minimizer_ref
             COMMAND ${CMAKE_BINARY_DIR}/bin/minimizer_ref --input ${CMAKE_SOURCE_DIR}/dataset/b64/test_64.nuc -o ${CMAKE_BINARY_DIR}/test_ref_output.txt
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    
    add_test(NAME test_minimizer_thr_v1
             COMMAND ${CMAKE_BINARY_DIR}/bin/minimizer_thr_v1 --input ${CMAKE_SOURCE_DIR}/dataset/b64/test_64.nuc -o ${CMAKE_BINARY_DIR}/test_thr_output.txt
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    
    add_test(NAME compare_outputs
             COMMAND ${CMAKE_COMMAND} -E compare_files
             ${CMAKE_BINARY_DIR}/test_ref_output.txt
             ${CMAKE_BINARY_DIR}/test_thr_output.txt)
endif()

message(STATUS "Configuration terminée. Pour compiler, exécutez 'make'")
message(STATUS "Pour installer: make install")
message(STATUS "Pour exécuter les tests: make test")
message(STATUS "Pour nettoyer: make clean")
