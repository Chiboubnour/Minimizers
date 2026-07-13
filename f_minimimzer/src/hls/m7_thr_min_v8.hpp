#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m7 : minimum glissant (template sur PAR_FACTOR), header-only.
//
//  Prologue : charge le warm-up de (WINDOW_SIZE-1) hashes dans "mem". Regime
//  permanent : win = (nouveaux, mem) ; pour chaque slot i, minimum des WINDOW_SIZE
//  hashes win[i .. i+WINDOW_SIZE-1] (extraction parallele). Masque valid relaye.
//
template<int PAR_FACTOR>
void thr_min_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  hash_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS       = 2 * SMER_SIZE;                          // un hash
    constexpr int HIST_BITS    = (WINDOW_SIZE - 1) * S_BITS;             // warm-up (bits)
    constexpr int WIN_BITS     = HIST_BITS + PAR_FACTOR * S_BITS;        // fenetre combinee
    constexpr int first_rounds = (WINDOW_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (WINDOW_SIZE - 1) % PAR_FACTOR;

    ap_uint<HIST_BITS> mem = 0;             // WINDOW_SIZE-1 hashes, plus ancien en bas

    ap_uint<S_BITS * PAR_FACTOR> in_word;
    ap_uint<         PAR_FACTOR> in_valid;

    //
    // Prologue : first_rounds paquets pleins + 1 paquet de n_last_round hashes.
    //
    for (int r = 0; r < first_rounds; r++) {
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        mem.range(S_BITS * PAR_FACTOR * (r + 1) - 1, S_BITS * PAR_FACTOR * r) = in_word;
    }
    in_word  = hash_stream_i.read();
    in_valid = hash_valid_i.read();
    mem.range(HIST_BITS - 1, S_BITS * PAR_FACTOR * first_rounds) = in_word.range(n_last_round * S_BITS - 1, 0);

    //
    // Regime permanent.
    //
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();

        if (in_valid == 0) {
            break;
        }

        const ap_uint<WIN_BITS> win = (in_word, mem);

        ap_uint<S_BITS * PAR_FACTOR> packed_out = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
#pragma HLS UNROLL
            ap_uint<S_BITS> m = win.range((i + 1) * S_BITS - 1, i * S_BITS);
            for (int w = 1; w < WINDOW_SIZE; w++) {
                const ap_uint<S_BITS> v = win.range((i + w + 1) * S_BITS - 1, (i + w) * S_BITS);
                if (v < m) m = v;
            }
            packed_out.range((i + 1) * S_BITS - 1, i * S_BITS) = m;
        }

        min_stream_o.write(packed_out);
        min_valid_o.write (in_valid);

        mem = win.range(WIN_BITS - 1, PAR_FACTOR * S_BITS);
    }

    min_stream_o.write(0);
    min_valid_o.write (0);
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
