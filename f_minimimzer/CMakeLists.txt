cmake_minimum_required(VERSION 3.21)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_BUILD_TYPE Release)

set (CMAKE_CXX_STANDARD_REQUIRED ON)

project(NB_SC_DECODER_EXPLORER LANGUAGES CXX)

option (SANITIZE_DECODERS "Add sanitization values" OFF)
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

add_executable(main
	./src/tools/compare.cpp
	./src/tools/TimeMeasure.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/tools/stats.cpp
	./src/main.cpp
	)
target_include_directories (main PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/common)

#######################################################################################################################################

project(test)
#
set(SOURCE_FILES
	./src/tools/compare.cpp
	./src/tools/TimeMeasure.cpp
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/tools/stats.cpp
	./src/test.cpp
)
#
add_executable(test ${SOURCE_FILES})

#######################################################################################################################################
# Testbench unitaire du module m1 (thr_reader), isole des autres tests.

add_executable(tb_m1
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/test_m1.cpp
	)

#######################################################################################################################################
# Testbench unitaire du module m2 (encodeur 2 bits), alimente par les sorties de m1.

add_executable(tb_m2
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/test_m2.cpp
	)

#######################################################################################################################################
# Testbench unitaire du module m3 (aligneur de s-mer), chaine m1 -> m2 -> m3.

add_executable(tb_m3
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/test_m3.cpp
	)

#######################################################################################################################################
# Testbench unitaire du module m4 (generateur de s-mers), chaine m1 -> m2 -> m3 -> m4.

add_executable(tb_m4
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/test_m4.cpp
	)

#######################################################################################################################################
# Testbench unitaire du module m5 (hash des s-mers), chaine m1 -> ... -> m5.

add_executable(tb_m5
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/hls/m5_thr_hash.cpp
	./src/test_m5.cpp
	)

#######################################################################################################################################
# Testbench unitaire du module m6 (aligneur de fenetre), chaine m1 -> ... -> m6.

add_executable(tb_m6
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/hls/m5_thr_hash.cpp
	./src/hls/m6_thr_adapter_smer.cpp
	./src/test_m6.cpp
	)

#######################################################################################################################################
# Testbench unitaire du module m7 (minimum glissant), chaine m1 -> ... -> m7.

add_executable(tb_m7
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/hls/m5_thr_hash.cpp
	./src/hls/m6_thr_adapter_smer.cpp
	./src/hls/m7_thr_min_v8.cpp
	./src/test_m7.cpp
	)

#######################################################################################################################################
# Testbench unitaire du module m8 (deduplication), chaine m1 -> ... -> m8.

add_executable(tb_m8
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/hls/m5_thr_hash.cpp
	./src/hls/m6_thr_adapter_smer.cpp
	./src/hls/m7_thr_min_v8.cpp
	./src/hls/m8_thr_dedup_v8.cpp
	./src/test_m8.cpp
	)

#######################################################################################################################################
# Testbench du module m9 (ecriture memoire finale), chaine complete m1 -> ... -> m9.

add_executable(tb_m9
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/hls/m5_thr_hash.cpp
	./src/hls/m6_thr_adapter_smer.cpp
	./src/hls/m7_thr_min_v8.cpp
	./src/hls/m8_thr_dedup_v8.cpp
	./src/hls/m9_thr_store_burst.cpp
	./src/test_m9.cpp
	)

#######################################################################################################################################
# Testbench top-level : le noyau minimizer() assemble (m1 -> ... -> m9) en dataflow.

add_executable(tb_minimizer
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/hls/m5_thr_hash.cpp
	./src/hls/m6_thr_adapter_smer.cpp
	./src/hls/m7_thr_min_v8.cpp
	./src/hls/m8_thr_dedup_v8.cpp
	./src/hls/m9_thr_store_burst.cpp
	./src/hls/minimizer.cpp
	./src/test_minimizer.cpp
	)

#######################################################################################################################################
# Tops de synthese concrets (minimizer_pf8/16/32/64) : compile-check des 4 instanciations.

add_library(minimizer_wrappers OBJECT
	./src/hls/minimizer_top.cpp
	)

#######################################################################################################################################
# Test de reference (golden vector) sur une sequence ADN concrete.

add_executable(tb_sequence
	./src/catch2v3/catch_amalgamated.cpp
	./src/hls/m1_thr_reader.cpp
	./src/hls/m2_thr_base_2b_conv.cpp
	./src/hls/m3_thr_base_adapter.cpp
	./src/hls/m4_thr_smer_gen.cpp
	./src/hls/m5_thr_hash.cpp
	./src/hls/m6_thr_adapter_smer.cpp
	./src/hls/m7_thr_min_v8.cpp
	./src/hls/m8_thr_dedup_v8.cpp
	./src/hls/m9_thr_store_burst.cpp
	./src/test_sequence.cpp
	)

#######################################################################################################################################
