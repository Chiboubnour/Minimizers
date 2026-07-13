#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m8 : deduplication des minimizers consecutifs identiques (template sur
//  PAR_FACTOR), header-only.
//
//  Pas de compaction : un doublon (== au dernier minimizer garde, etat "last"
//  conserve entre paquets) a juste son bit valid mis a 0, a sa position. Un paquet
//  entierement masque n'est pas emis.
//
template<int PAR_FACTOR>
void thr_dedup_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  dedup_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  dedup_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * SMER_SIZE;   // un minimizer

    ap_uint<S_BITS> last     = 0;
    bool            has_last = false;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<2 * SMER_SIZE * PAR_FACTOR> in_word  = min_stream_i.read();
        const ap_uint<                PAR_FACTOR> in_valid = min_valid_i.read();

        if (in_valid == 0) {
            break;
        }

        ap_uint<PAR_FACTOR> out_valid = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (in_valid[i]) {
                const ap_uint<S_BITS> v = in_word.range((i + 1) * S_BITS - 1, i * S_BITS);
                if (!has_last || (v != last)) {
                    out_valid[i] = 1;       // conserve a sa position
                    last     = v;
                    has_last = true;
                }
            }
        }

        if (out_valid != 0) {
            dedup_stream_o.write(in_word);      // donnees inchangees
            dedup_valid_o.write (out_valid);    // masque creux
        }
    }

    dedup_stream_o.write(0);
    dedup_valid_o.write (0);
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
