cmake_minimum_required(VERSION 3.21)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_BUILD_TYPE Release)

set (CMAKE_CXX_STANDARD_REQUIRED ON)

project(NB_SC_DECODER_EXPLORER LANGUAGES CXX)

option (SANITIZE_DECODERS "Add sanitization values" ON)
if (SANITIZE_DECODERS)
	set (SANITIZE_CXXFLAGS "-fsanitize=address;-fsanitize=undefined;-fno-omit-frame-pointer;-ggdb")
	set (SANITIZE_LDFLAGS  "-fsanitize=address;-fsanitize=undefined;-fno-omit-frame-pointer;-ggdb")
endif ()

add_compile_options (-march=native -mtune=native -O3 -g0 ${SANITIZE_CXXFLAGS})
add_compile_options (-Wall -Wextra)
add_compile_options (-I../include)
# to remove warnings from ap_int
add_compile_options (-Wno-unused-parameter -Wno-deprecated-copy -Wno-int-in-bool-context -Wno-type-limits -Wno-strict-aliasing -Wno-extra -Wno-unused-label -Wno-unknown-pragmas)
add_link_options (${SANITIZE_LDFLAGS})

#######################################################################################################################################

add_executable(hls_test
		./hls/fasta_parser.cpp
		./src/main.cpp
	)
target_include_directories (hls_test PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/common)

#######################################################################################################################################
