#include <ap_int.h>
#include <hls_stream.h>

void adapter_hls(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<64>>& out_stream,
    hls::stream<ap_uint<8>>&  out_valid
){
#pragma HLS INLINE off

    ap_uint<128> buffer = 0;
    ap_uint<8> valid_count = 0;


    while (true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> in_word  = base_stream.read();
        ap_uint<8>  in_valid = base_valid.read();
        if (in_valid == 0){

            if(valid_count > 0){
                ap_uint<64> out_word = buffer.range(63,0);
                out_stream.write(out_word);
                out_valid.write(valid_count);
            }


            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        buffer.range(valid_count*8 + 63, valid_count*8) = in_word;
        valid_count += in_valid;

        if(valid_count >= 8){

            ap_uint<64> out_word = buffer.range(63,0);
            out_stream.write(out_word);
            out_valid.write(8);

            buffer >>= 64;
            valid_count -= 8;
        }
    }
}

