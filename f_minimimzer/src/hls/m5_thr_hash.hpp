#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m5 : hachage des s-mers canoniques (template sur PAR_FACTOR), header-only.
//
//  Pour chaque s-mer (2*SMER_SIZE bits), son hash sur le meme nombre de bits
//  (melange 64 bits de Wang tronque). Masque valid relaye inchange.
//
template <int resu_size>
inline ap_uint<resu_size> hash_u64(ap_uint<64> key) {
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return key;
}
//
//
template<int PAR_FACTOR>
void thr_hash(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  hash_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * SMER_SIZE;   // largeur d'un s-mer (et de son hash)

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_smer = smer_stream_i.read();
        const ap_uint<                PAR_FACTOR> valid       = smer_valid_i.read();

        if (valid == 0) {
            hash_stream_o.write(0);
            hash_valid_o.write (0);
            break;
        }

        ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_hash = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
#pragma HLS UNROLL
            const ap_uint<S_BITS> smer = packed_smer.range((i + 1) * S_BITS - 1, i * S_BITS);
            const ap_uint<S_BITS> h    = hash_u64<S_BITS>((ap_uint<64>)smer);
            packed_hash.range((i + 1) * S_BITS - 1, i * S_BITS) = h;
        }

        hash_stream_o.write(packed_hash);
        hash_valid_o.write (valid);
    }
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
