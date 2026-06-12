#pragma once
#include "./header.hpp"
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
extern void thr_smer_gen(
    hls::stream<ap_uint<2 *             PAR_FACTOR>>& base_stream_i,
    hls::stream<ap_uint<                PAR_FACTOR>>& base_valid_i,
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_valid_o
);
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//