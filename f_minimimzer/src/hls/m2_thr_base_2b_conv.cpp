#include "./m2_thr_base_2b_conv.hpp"
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
void thr_adapter_hls(
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

        const ap_uint<64> in_word  = base_stream_i.read();
        const ap_uint<8>  in_valid = base_valid_i.read();

        for(int i = 0; i < PAR_FACTOR; i += 1)
        {
            const ap_uint<8> nucl  = in_word.range(8 * i + 7, 8 * i);
            const ap_uint<2> enc2b = (nucl >> 1) & 3;
            buffer(2 * (i+1) - 1, 2 * i) = enc2b; 
        }

        base_stream_o.write( buffer   );
        base_valid_o.write ( in_valid );

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