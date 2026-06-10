//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
template<int S_BASE>
void hash_thread(
    hls::stream<ap_uint<Parallel_factor*2*S_BASE>>& smer_stream,
    hls::stream<ap_uint<8>>&  smer_valid,
    hls::stream<ap_uint<Parallel_factor*2*S_BASE>>& hash_stream,
    hls::stream<ap_uint<8>>&  hash_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<Parallel_factor*S_BITS> packed_smer = smer_stream.read();
        ap_uint<8> valid = smer_valid.read();

        if(valid == 0){
            hash_stream.write(0);
            hash_valid.write(0);
            break;
        }

        ap_uint<Parallel_factor*S_BITS> packed_hash = 0;

    HASH_LOOP:
        for(int i = 0; i < Parallel_factor; i++){

            if(i < valid){
                ap_uint<S_BITS> smer = packed_smer.range((i+1)*S_BITS-1, i*S_BITS);
                ap_uint<S_BITS> h    = hash_u64<S_BITS>((ap_uint<64>)smer);
                packed_hash.range((i+1)*S_BITS-1, i*S_BITS) = h;
            }
        }

        hash_stream.write(packed_hash);
        hash_valid.write(valid);
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