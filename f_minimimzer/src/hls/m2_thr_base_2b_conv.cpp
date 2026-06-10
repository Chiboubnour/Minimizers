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
static void thr_adapter_hls(
    hls::stream<ap_uint<64>>& base_stream_i,
    hls::stream<ap_uint<8>>&  base_valid_i,
    hls::stream<ap_uint<16>>& base_stream_o,
    hls::stream<ap_uint<8>>&  base_valid_o
){
#pragma HLS INLINE off

    // buffer 128 bits pour absorber 2 mots de 64 bits (64 bases 2 bits)
    ap_uint<16> buffer;
    ap_uint<8>  valid_count;

    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> in_word  = base_stream_i.read();
        ap_uint<8>  in_valid = base_valid_i.read();

        for(int i = 0; i < PAR_FACTOR; i += 1)
        {
            buffer = 
        }

        out_stream_o.write( buffer );
        out_valid_o.write( in_valid );

        if (in_valid == 0)
        {
            break;
        }
    }
}
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//