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
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<64>>& out_stream,
    hls::stream<ap_uint<8>>&  out_valid
){
#pragma HLS INLINE off

    // buffer 128 bits pour absorber 2 mots de 64 bits (64 bases 2 bits)
    ap_uint<128> buffer      = 0;
    ap_uint<8>   valid_count = 0;

    while (true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> in_word  = base_stream.read();
        ap_uint<8>  in_valid = base_valid.read();

        if (in_valid == 0){

            if(valid_count > 0){
                out_stream.write(buffer.range(63, 0));
                out_valid.write(valid_count);
            }

            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        // chaque base occupe 2 bits : on décale de valid_count*2 bits
        buffer.range(valid_count*2 + 63, valid_count*2) = in_word;
        valid_count += in_valid;

        if(valid_count >= Parallel_factor){

            out_stream.write(buffer.range(63, 0));
            out_valid.write(Parallel_factor);

            buffer >>= 64;
            valid_count -= Parallel_factor;
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