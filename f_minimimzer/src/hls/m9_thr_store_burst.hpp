#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m9 : ecriture memoire des minimizers (template sur PAR_FACTOR), header-only.
//
//  Lit le flux a masque creux de m8 ; n'ecrit en memoire que les minimizers
//  valides, tasses densement (adresses consecutives = burst). nElem = leur nombre.
//
template<int PAR_FACTOR>
void thr_store_burst(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  minz_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  minz_valid_i,
    ap_uint<64>*                                        tab_hash,
    ap_uint<64>*                                        nElem
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * SMER_SIZE;   // un minimizer

    ap_uint<64> cnt = 0;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed = minz_stream_i.read();
        const ap_uint<                PAR_FACTOR> valid  = minz_valid_i.read();

        if (valid == 0) {
            break;                          // paquet terminal
        }

        // collecte des minimizers valides du paquet (masque creux -> dense)
        ap_uint<S_BITS> burst[PAR_FACTOR];
#pragma HLS ARRAY_PARTITION variable=burst complete
        int n = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (valid[i]) {
                burst[n] = packed.range((i + 1) * S_BITS - 1, i * S_BITS);
                n++;
            }
        }

        // ecriture burst (adresses consecutives)
        for (int j = 0; j < n; j++) {
#pragma HLS PIPELINE II=1
            tab_hash[cnt + j] = burst[j];   // 2*SMER_SIZE bits -> zero-extension sur 64 bits
        }
        cnt += n;
    }

    *nElem = cnt;
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
