#include <ap_int.h>
#include <hls_stream.h>

#define S_BASE 5

template <int resu_size>
inline ap_uint<resu_size> hash_u64(ap_uint<64> key) {
#pragma HLS INLINE
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return key;
}

void smer_thread_hls(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<16*S_BASE>>& out_stream,
    hls::stream<ap_uint<8>>&  out_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;
    const ap_uint<S_BITS> MASK = ap_uint<S_BITS>(-1);

    ap_uint<S_BITS> fwd_state = 0;
    ap_uint<S_BITS> rev_state = 0;
    ap_uint<32> base_cnt = 0;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1

        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if(valid == 0){
            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        ap_uint<8*S_BITS> packed_out = 0;

        ap_uint<S_BITS> fwd_local = fwd_state;
        ap_uint<S_BITS> rev_local = rev_state;

        int out_idx = 0;

        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            ap_uint<8> nucl = word.range(8*i+7,8*i);
            ap_uint<2> c = (nucl >> 1) & 3;
            ap_uint<2> comp = 0x2 ^ c;

            fwd_local <<= 2;
            fwd_local(1,0) = c;

            rev_local >>= 2;
            rev_local(S_BITS-1,S_BITS-2) = comp;

            base_cnt++;

            if(base_cnt >= S_BASE){
                ap_uint<S_BITS> canon =
                    (fwd_local < rev_local) ? fwd_local : rev_local;

                ap_uint<S_BITS> h =
                    hash_u64<S_BITS>((ap_uint<64>)canon);

                packed_out.range((out_idx+1)*S_BITS-1, out_idx*S_BITS) = h;
                out_idx++;
            }
        }

        fwd_state = fwd_local & MASK;
        rev_state = rev_local;

        if(out_idx > 0){
            out_stream.write(packed_out);
            out_valid.write(out_idx);
        }
    }
}

void smer_generate_thread(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<16*S_BASE>>& smer_stream,
    hls::stream<ap_uint<8>>&  smer_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;
    const ap_uint<S_BITS> MASK = ap_uint<S_BITS>(-1);

    ap_uint<S_BITS> fwd_state = 0;
    ap_uint<S_BITS> rev_state = 0;
    ap_uint<32> base_cnt = 0;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1

        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if(valid == 0){
            smer_stream.write(0);
            smer_valid.write(0);
            break;
        }

        ap_uint<8*S_BITS> packed_out = 0;

        ap_uint<S_BITS> fwd_local = fwd_state;
        ap_uint<S_BITS> rev_local = rev_state;

        int out_idx = 0;

        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            ap_uint<8> nucl = word.range(8*i+7,8*i);
            ap_uint<2> c = (nucl >> 1) & 3;
            ap_uint<2> comp = 0x2 ^ c;

            fwd_local <<= 2;
            fwd_local(1,0) = c;

            rev_local >>= 2;
            rev_local(S_BITS-1,S_BITS-2) = comp;

            base_cnt++;

            if(base_cnt >= S_BASE){
                ap_uint<S_BITS> canon =
                    (fwd_local < rev_local) ? fwd_local : rev_local;

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

void hash_thread(
    hls::stream<ap_uint<16*S_BASE>>& smer_stream,
    hls::stream<ap_uint<8>>&  smer_valid,
    hls::stream<ap_uint<16*S_BASE>>& hash_stream,
    hls::stream<ap_uint<8>>&  hash_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1

        ap_uint<8*S_BITS> packed_smer = smer_stream.read();
        ap_uint<8> valid = smer_valid.read();

        if(valid == 0){
            hash_stream.write(0);
            hash_valid.write(0);
            break;
        }

        ap_uint<8*S_BITS> packed_hash = 0;

        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            if(i < valid){
                ap_uint<S_BITS> smer =
                    packed_smer.range((i+1)*S_BITS-1, i*S_BITS);

                ap_uint<S_BITS> h =
                    hash_u64<S_BITS>((ap_uint<64>)smer);

                packed_hash.range((i+1)*S_BITS-1, i*S_BITS) = h;
            }
        }

        hash_stream.write(packed_hash);
        hash_valid.write(valid);
    }
}
