#pragma once
#include "./header.hpp"
//
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
ap_uint<8> valid_bits(const int reste)
{
    ap_uint<8> res;
    switch (reste)
    {
        case 0:  res = 0x00; break;
        case 1:  res = 0x01; break;
        case 2:  res = 0x03; break;
        case 3:  res = 0x07; break;
        case 4:  res = 0x0F; break;
        case 5:  res = 0x1F; break;
        case 6:  res = 0x3F; break;
        case 7:  res = 0x7F; break;
        default: res = 0xFF; break;
    }
    return res;
}

//
//
static void thr_reader(
    const ap_uint< MEM_WIDTH >* base_ptr_i,
    const ap_uint<        64 >  n_bases_i,
    hls::stream< ap_uint<64> >& base_stream_o,
    hls::stream< ap_uint< 8> >& base_valid_o
) {
#pragma HLS INLINE off
    //
    // chaque mot 64 bits contient 32 bases encodées sur 2 bits
    //
    ap_uint<64> n_words  = (n_bases_i + PAR_FACTOR - 1) / PAR_FACTOR;
    ap_uint<64> base_cnt =  n_bases_i;
    for (ap_uint<64> w = 0; w < n_words; w++)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<64> word = base_ptr_i[w];
        ap_uint<8> valid = valid_bits( base_cnt );

        base_stream_o.write(word);
        base_valid_o.write(valid);

        base_cnt -= PAR_FACTOR;
    }
    base_stream.write(0);
    base_valid.write(0);
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