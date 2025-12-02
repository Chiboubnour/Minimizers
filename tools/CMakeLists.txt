cmake_minimum_required(VERSION 3.15)
project(Minimizers)

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_BUILD_TYPE  Release)
#SET(CMAKE_BUILD_TYPE  Debug)

option(BUILD_LINUX_INTEL  "Build on Linux with INTEL processor" OFF) # JetSon Orin Nano (ENSSAT)
option(BUILD_MACOS_ARM    "Build on macOs with ARM.  processor" ON)  # MACOS INTEL-Mx


if(BUILD_LINUX_INTEL)
    add_compile_options(-Wall -Wextra -Werror)
    add_compile_options(-O3 -march=native -mtune=native)
    #add_compile_options(-fopenmp)
    #add_link_options(-fopenmp)
endif()

if(BUILD_MACOS_ARM)
    # TO HAVE A NATIVE OPENMP SUPPORT
    #SET(CMAKE_C_COMPILER   /opt/homebrew/opt/llvm/bin/clang)
    #SET(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)
    add_compile_options(-Wall -Wextra)
    add_compile_options(-O3 -mcpu=native -mtune=native)
    #add_compile_options(-fopenmp)
    #add_link_options(-fopenmp)
endif()


# ---- EXEC MONO ----
add_executable(nucl_extractor
    ../apps/nucl_extractor.cpp
    ../src/files/stream_reader_library.cpp
    ../src/files/gz/reader/stream_gz_reader.cpp
    ../src/files/gz/writer/stream_gz_writer.cpp
    ../src/files/stream_writer_library.cpp
    ../src/files/lz4/reader/stream_lz4_reader.cpp
    ../src/files/lz4/writer/stream_lz4_writer.cpp
    ../src/files/bz2/reader/stream_bz2_reader.cpp
    ../src/files/bz2/writer/stream_bz2_writer.cpp
    ../src/files/raw/reader/stream_raw_reader.cpp
    ../src/files/raw/writer/stream_raw_writer.cpp
    ../src/stats/file_stats.cpp
    # File parsers
    ../src/parsers/parser_fasta.cpp
    ../src/parsers/parser_fastq.cpp
    # LZ4 lib
    ../src/libs/lz4/lz4hc.c
    ../src/libs/lz4/lz4frame.c
    ../src/libs/lz4/xxhash.c
    ../src/libs/lz4/lz4.c
    ../src/libs/lz4/lz4file.c    
)
target_link_libraries(nucl_extractor z)
target_link_libraries(nucl_extractor bz2)

# ---- EXEC MONO ----
add_executable(nucl_encode
    ../apps/nucl_encode.cpp
    ../src/files/stream_reader_library.cpp
    ../src/files/gz/reader/stream_gz_reader.cpp
    ../src/files/gz/writer/stream_gz_writer.cpp
    ../src/files/stream_writer_library.cpp
    ../src/files/lz4/reader/stream_lz4_reader.cpp
    ../src/files/lz4/writer/stream_lz4_writer.cpp
    ../src/files/bz2/reader/stream_bz2_reader.cpp
    ../src/files/bz2/writer/stream_bz2_writer.cpp
    ../src/files/raw/reader/stream_raw_reader.cpp
    ../src/files/raw/writer/stream_raw_writer.cpp
    ../src/stats/file_stats.cpp
    # LZ4 lib
    ../src/libs/lz4/lz4hc.c
    ../src/libs/lz4/lz4frame.c
    ../src/libs/lz4/xxhash.c
    ../src/libs/lz4/lz4.c
    ../src/libs/lz4/lz4file.c    
)
target_link_libraries(nucl_encode z)
target_link_libraries(nucl_encode bz2)

#target_link_libraries(minimizer_mono PRIVATE common_interface)
