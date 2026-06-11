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
extern void m3_thr_base_adapter(
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& base_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid,
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& out_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& out_valid
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