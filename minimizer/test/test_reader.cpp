/*
 *  Copyright (c) 2022 Bertrand LE GAL
 *
 *  This software is provided 'as-is', without any express or
 *  implied warranty. In no event will the authors be held
 *  liable for any damages arising from the use of this software.
 *
 *  Permission is granted to anyone to use this software for any purpose,
 *  including commercial applications, and to alter it and redistribute
 *  it freely, subject to the following restrictions:
 *
 *  1. The origin of this software must not be misrepresented;
 *  you must not claim that you wrote the original software.
 *  If you use this software in a product, an acknowledgment
 *  in the product documentation would be appreciated but
 *  is not required.
 *
 *  2. Altered source versions must be plainly marked as such,
 *  and must not be misrepresented as being the original software.
 *
 *  3. This notice may not be removed or altered from any
 *  source distribution.
 *
 */
//
//
#define CATCH_CONFIG_MAIN
#include "./catch2v3/catch_amalgamated.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <ap_int.h>
#include "functions.hpp"

static void reader_hls(
    const uint64_t* packed_sequence,
    ap_uint<64> n_bases,
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid
) {
#pragma HLS INLINE off

    ap_uint<64> n_words = (n_bases + 7) >> 3;

    for (ap_uint<64> w=0; w<n_words; w++) {

        ap_uint<64> word = packed_sequence[w];
        ap_uint<8> valid;
        if(w == n_words-1 && (n_bases & 7)) {
            valid = n_bases & 7; 
        } else {
            valid = 8;
        }

        base_stream.write(word);
        base_valid.write(valid);
    }

    base_stream.write(0);
    base_valid.write(0);
}

/* ============================================================
   Helper : vérifie qu’un stream est vide après lecture
============================================================ */
static void check_stream_empty(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>& valid_stream)
{
    REQUIRE(base_stream.empty());
    REQUIRE(valid_stream.empty());
}

/* ============================================================
   R1 — Séquence multiple de 8
============================================================ */
TEST_CASE("R1_reader_nominal_64_bases", "[reader_hls]")
{
    const uint64_t seq[] = {
        0x5441414741414154,
        0x4743544154474754,
        0x5447474754434341,
        0x4747544341544341,
        0x5441414741414154,
        0x4743544154474754,
        0x5447474754434341,
        0x4747544341544341
    };

    hls::stream<ap_uint<64>> base_o;
    hls::stream<ap_uint<8>>  valid_o;

    reader_hls(seq, 64, base_o, valid_o);

    for(int i=0;i<8;i++)
        REQUIRE(base_o.read() == seq[i]);

    REQUIRE(base_o.read() == 0);

    for(int i=0;i<8;i++)
        REQUIRE(valid_o.read() == 8);

    REQUIRE(valid_o.read() == 0);

    check_stream_empty(base_o, valid_o);
}

/* ============================================================
   R2 — Mot final partiel
============================================================ */
TEST_CASE("R2_reader_partial_last_word_66_bases", "[reader_hls]")
{
    const uint64_t seq[] = {
        0x5441414741414154,
        0x4743544154474754,
        0x5447474754434341,
        0x4747544341544341,
        0x5441414741414154,
        0x4743544154474754,
        0x5447474754434341,
        0x4747544341544341,
        0x0000000000004341
    };

    hls::stream<ap_uint<64>> base_o;
    hls::stream<ap_uint<8>>  valid_o;

    reader_hls(seq, 66, base_o, valid_o);

    for(int i=0;i<8;i++)
        REQUIRE(base_o.read() == seq[i]);

    REQUIRE(base_o.read() == seq[8]);
    REQUIRE(base_o.read() == 0);

    for(int i=0;i<8;i++)
        REQUIRE(valid_o.read() == 8);

    REQUIRE(valid_o.read() == 2);
    REQUIRE(valid_o.read() == 0);

    check_stream_empty(base_o, valid_o);
}

/* ============================================================
   R3 — Séquence vide
============================================================ */
TEST_CASE("R3_reader_empty_sequence", "[reader_hls]")
{
    uint64_t dummy = 0;

    hls::stream<ap_uint<64>> base_o;
    hls::stream<ap_uint<8>>  valid_o;

    reader_hls(&dummy, 0, base_o, valid_o);

    REQUIRE(base_o.read() == 0);
    REQUIRE(valid_o.read() == 0);

    check_stream_empty(base_o, valid_o);
}

/* ============================================================
   R4 — Une seule base
============================================================ */
TEST_CASE("R4_reader_single_base", "[reader_hls]")
{
    const uint64_t seq[] = { 0x0000000000000041 };

    hls::stream<ap_uint<64>> base_o;
    hls::stream<ap_uint<8>>  valid_o;

    reader_hls(seq, 1, base_o, valid_o);

    REQUIRE(base_o.read() == seq[0]);
    REQUIRE(base_o.read() == 0);

    REQUIRE(valid_o.read() == 1);
    REQUIRE(valid_o.read() == 0);

    check_stream_empty(base_o, valid_o);
}

/* ============================================================
   R5 — Taille SMER-1 
============================================================ */
TEST_CASE("R5_reader_boundary_smer_minus_1", "[reader_hls]")
{
    const uint64_t seq[] = { 0x4141414141414141 };

    hls::stream<ap_uint<64>> base_o;
    hls::stream<ap_uint<8>>  valid_o;

    reader_hls(seq, 7, base_o, valid_o);

    REQUIRE(base_o.read() == seq[0]);
    REQUIRE(base_o.read() == 0);

    REQUIRE(valid_o.read() == 7);
    REQUIRE(valid_o.read() == 0);

    check_stream_empty(base_o, valid_o);
}

/* ============================================================
   R6 — Boundary SMER
============================================================ */
TEST_CASE("R6_reader_boundary_smer_exact", "[reader_hls]")
{
    const uint64_t seq[] = {
        0xAAAAAAAAAAAAAAAA,
        0xBBBBBBBBBBBBBBBB,
        0xCCCCCCCCCCCCCCCC
    };

    hls::stream<ap_uint<64>> base_o;
    hls::stream<ap_uint<8>>  valid_o;

    reader_hls(seq, 24, base_o, valid_o);

    for(int i=0;i<3;i++)
        REQUIRE(base_o.read() == seq[i]);

    REQUIRE(base_o.read() == 0);

    for(int i=0;i<3;i++)
        REQUIRE(valid_o.read() == 8);

    REQUIRE(valid_o.read() == 0);

    check_stream_empty(base_o, valid_o);
}

/* ============================================================
   R7 — Longue séquence (robustesse index)
============================================================ */
TEST_CASE("R7_reader_long_sequence_stress", "[reader_hls]")
{
    constexpr int NWORDS = 128;
    uint64_t seq[NWORDS];

    for(int i=0;i<NWORDS;i++)
        seq[i] = i;

    hls::stream<ap_uint<64>> base_o;
    hls::stream<ap_uint<8>>  valid_o;

    reader_hls(seq, NWORDS*8, base_o, valid_o);

    for(int i=0;i<NWORDS;i++)
        REQUIRE(base_o.read() == seq[i]);

    REQUIRE(base_o.read() == 0);

    for(int i=0;i<NWORDS;i++)
        REQUIRE(valid_o.read() == 8);

    REQUIRE(valid_o.read() == 0);

    check_stream_empty(base_o, valid_o);
}
