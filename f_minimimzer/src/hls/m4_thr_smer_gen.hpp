#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m4 : generateur de s-mers canoniques (template sur PAR_FACTOR), header-only.
//
//  Convention "plus RECENTE en bas" : dans un s-mer, la base la plus recente est
//  dans les bits de poids faible (necessaire pour reproduire le hash de reference).
//
//  Prologue : charge le warm-up de (SMER_SIZE-1) bases dans "memory" (recent en
//  bas). Regime permanent : on reordonne les nouvelles bases (recente en bas),
//  on forme win = (memory, new_rev), et on extrait PAR_FACTOR s-mers par plages
//  (sans dependance entre slots), puis min(forward, reverse-complement).
//
template<int PAR_FACTOR>
void thr_smer_gen(
    hls::stream<ap_uint<2 *             PAR_FACTOR>>& base_stream_i,
    hls::stream<ap_uint<                PAR_FACTOR>>& base_valid_i,
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_valid_o
){
#pragma HLS INLINE off

    constexpr int SMER_BITS    = 2 * SMER_SIZE;                          // 2*19 : un s-mer
    constexpr int HIST_BITS    = 2 * (SMER_SIZE - 1);                    // warm-up (bits)
    constexpr int WIN_BITS     = HIST_BITS + 2 * PAR_FACTOR;             // fenetre combinee
    constexpr int first_rounds = (SMER_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (SMER_SIZE - 1) % PAR_FACTOR;

    ap_uint<HIST_BITS> memory = 0;          // SMER_SIZE-1 bases, plus recente en bas

    ap_uint<2 * PAR_FACTOR> in_word;
    ap_uint<    PAR_FACTOR> in_valid;

    //
    // Prologue : charge les (SMER_SIZE-1) bases de warm-up, base par base (recente
    // en bas). first_rounds paquets pleins + 1 paquet de n_last_round bases.
    //
    for (int r = 0; r < first_rounds; r++) {
        in_word  = base_stream_i.read();
        in_valid = base_valid_i.read();
        for (int i = 0; i < PAR_FACTOR; i++) {
            memory = (memory << 2);
            memory(1, 0) = in_word.range(2 * i + 1, 2 * i);
        }
    }
    in_word  = base_stream_i.read();
    in_valid = base_valid_i.read();
    for (int i = 0; i < n_last_round; i++) {
        memory = (memory << 2);
        memory(1, 0) = in_word.range(2 * i + 1, 2 * i);
    }

    //
    // Regime permanent.
    //
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        in_word  = base_stream_i.read();
        in_valid = base_valid_i.read();

        if (in_valid == 0) {
            break;
        }
        //
        // reordonne les PAR_FACTOR nouvelles bases : "ancienne en bas" (in_word)
        // -> "recente en bas" (new_rev). Cablage, sans dependance portee.
        //
        ap_uint<2 * PAR_FACTOR> new_rev = 0;
        for (int j = 0; j < PAR_FACTOR; j++) {
#pragma HLS UNROLL
            new_rev.range(2 * (PAR_FACTOR - 1 - j) + 1, 2 * (PAR_FACTOR - 1 - j)) = in_word.range(2 * j + 1, 2 * j);
        }
        //
        // fenetre combinee : historique (haut) + nouvelles bases (bas), recente en bas
        //
        const ap_uint<WIN_BITS> win = (memory, new_rev);

        ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_out = 0;
        for (int i = 0; i < PAR_FACTOR; i++)
        {
#pragma HLS UNROLL
            // s-mer (forward) du slot i : recente en bas
            const ap_uint<SMER_BITS> fwd = win.range(2 * (PAR_FACTOR - 1 - i) + 2 * SMER_SIZE - 1,
                                                     2 * (PAR_FACTOR - 1 - i));
            // reverse-complement : ordre inverse + complement (A<->T, C<->G => XOR 0b10)
            ap_uint<SMER_BITS> rev = 0;
            for (int t = 0; t < SMER_SIZE; t++)
            {
#pragma HLS UNROLL
                const ap_uint<2> base = fwd.range(2 * t + 1, 2 * t);
                rev.range(2 * (SMER_SIZE - 1 - t) + 1, 2 * (SMER_SIZE - 1 - t)) = base ^ ap_uint<2>(0x2);
            }
            const ap_uint<SMER_BITS> canon = (fwd < rev) ? fwd : rev;
            packed_out.range((i + 1) * SMER_BITS - 1, i * SMER_BITS) = canon;
        }

        smer_stream_o.write(packed_out);
        smer_valid_o.write (in_valid);
        //
        // memorise les (SMER_SIZE-1) bases les plus recentes (le bas de win)
        //
        memory = win.range(HIST_BITS - 1, 0);
    }

    smer_stream_o.write(0);
    smer_valid_o.write (0);
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
