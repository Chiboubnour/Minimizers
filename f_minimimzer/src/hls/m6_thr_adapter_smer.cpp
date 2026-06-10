//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
template<int S_BITS>
static void thr_adapter_smer(
    hls::stream<ap_uint<Parallel_factor*S_BITS>>& base_stream,
    hls::stream<ap_uint<8>>&                      base_valid,
    hls::stream<ap_uint<Parallel_factor*S_BITS>>& out_stream,
    hls::stream<ap_uint<8>>&                      out_valid
)
{
#pragma HLS INLINE off

    ap_uint<2*Parallel_factor*S_BITS> buffer      = 0;
    ap_uint<8>                        valid_count  = 0;

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto in_word  = base_stream.read();
        auto in_valid = base_valid.read();

        if(in_valid == 0){

            if(valid_count > 0){
                out_stream.write(buffer.range(Parallel_factor*S_BITS-1, 0));
                out_valid.write(valid_count);
            }

            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        buffer |= ((ap_uint<2*Parallel_factor*S_BITS>)in_word) << (valid_count*S_BITS);
        valid_count += in_valid;

        if(valid_count >= Parallel_factor){

            out_stream.write(buffer.range(Parallel_factor*S_BITS-1, 0));
            out_valid.write(Parallel_factor);

            buffer >>= Parallel_factor*S_BITS;
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