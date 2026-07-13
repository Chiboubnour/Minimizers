#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m2 : encodeur 2 bits (template sur PAR_FACTOR), header-only.
//
//  Convertit PAR_FACTOR bases ASCII (mot de 8*PAR_FACTOR bits) en PAR_FACTOR codes
//  de 2 bits (mot de 2*PAR_FACTOR bits), via (nucl >> 1) & 3 : A=00, C=01, T=10,
//  G=11. Le masque valid est relaye inchange. Le paquet terminal est propage.
//
template<int PAR_FACTOR>
void thr_adapter_hls(
    hls::stream<ap_uint< 8 * PAR_FACTOR >>& base_stream_i,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid_i,
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& base_stream_o,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid_o
){
#pragma HLS INLINE off

    ap_uint<2 * PAR_FACTOR> buffer = 0;

    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<8 * PAR_FACTOR> in_word  = base_stream_i.read();
        const ap_uint<    PAR_FACTOR> in_valid = base_valid_i.read();

        for (int i = 0; i < PAR_FACTOR; i++)
        {
            const ap_uint<8> nucl  = in_word.range(8 * i + 7, 8 * i);
            const ap_uint<2> enc2b = (nucl >> 1) & 3;
            buffer.range(2 * i + 1, 2 * i) = enc2b;
        }

        base_stream_o.write(buffer);
        base_valid_o.write (in_valid);

        if (in_valid == 0)
        {
            break;
        }
    }
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
