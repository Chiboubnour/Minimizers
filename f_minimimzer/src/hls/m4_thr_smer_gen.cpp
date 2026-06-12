#include "./m4_thr_smer_gen.hpp"
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
void thr_smer_gen(
    hls::stream<ap_uint<2 *             PAR_FACTOR>>& base_stream_i,
    hls::stream<ap_uint<                PAR_FACTOR>>& base_valid_i,
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_valid_o
){
#pragma HLS INLINE off

    constexpr int SMER_BITS = 2 * S_BASE;
    const ap_uint<SMER_BITS> MASK = ap_uint<S_BITS>(-1);

    ap_uint<S_BITS> fwd_state = 0;
    ap_uint<S_BITS> rev_state = 0;
    ap_uint<32>     base_cnt  = 0;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if(valid == 0){
            smer_stream.write(0);
            smer_valid.write(0);
            break;
        }

        ap_uint<Parallel_factor*S_BITS> packed_out = 0;
        ap_uint<S_BITS> fwd_local = fwd_state;
        ap_uint<S_BITS> rev_local = rev_state;
        int out_idx = 0;

    BASE_LOOP:
        for(int i = 0; i < Parallel_factor; i++){

            // lecture directe des 2 bits de la base i dans le mot 64 bits
            ap_uint<2> c    = word.range(2*i+1, 2*i);
            ap_uint<2> comp = 0x2 ^ c;

            fwd_local <<= 2;
            fwd_local(1,0) = c;

            rev_local >>= 2;
            rev_local(S_BITS-1, S_BITS-2) = comp;

            base_cnt++;

            if(base_cnt >= S_BASE){
                ap_uint<S_BITS> canon = (fwd_local < rev_local) ? fwd_local : rev_local;
                packed_out.range((out_idx+1)*S_BITS-1, out_idx*S_BITS) = canon;
                out_idx++;
            }
        }

        fwd_state = fwd_local & MASK;
        rev_state = rev_local;

        if(out_idx > 0){
            smer_stream.write(packed_out);
            smer_valid.write(out_idx);
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
