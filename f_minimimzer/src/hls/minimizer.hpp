#pragma once
#include "./m1_thr_reader.hpp"
#include "./m2_thr_base_2b_conv.hpp"
#include "./m3_thr_base_adapter.hpp"
#include "./m4_thr_smer_gen.hpp"
#include "./m5_thr_hash.hpp"
#include "./m6_thr_adapter_smer.hpp"
#include "./m7_thr_min_v8.hpp"
#include "./m8_thr_dedup_v8.hpp"
#include "./m9_thr_store_burst.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  Top-level : calcul des minimizers d'une sequence ADN, en dataflow.
//  Template sur PAR_FACTOR : instancier minimizer<N> pour generer un circuit avec
//  un facteur de parallelisme N (8, 16, 32, 64...).
//
//  m1 lecture -> m2 encodage 2 bits -> m3 align s-mer -> m4 s-mers canoniques
//  -> m5 hash -> m6 align fenetre -> m7 min glissant -> m8 dedup -> m9 store.
//
//    packed_sequence : bases ASCII packees (PAR_FACTOR bases / mot de 8*PAR_FACTOR bits)
//    n_bases         : nombre de bases a traiter
//    tab_hash        : minimizers distincts (1 par mot de 64 bits)
//    nMinizrs        : nombre de minimizers ecrits
//
template<int PAR_FACTOR>
void minimizer(
    const ap_uint< 8 * PAR_FACTOR >* packed_sequence,
    ap_uint<64>*                     tab_hash,
    ap_uint<64>*                     nMinizrs,
    ap_uint<64>                      n_bases
){
    // Les pragmas d'INTERFACE (m_axi / s_axilite) sont portes par les wrappers
    // concrets (cf. minimizer_top.hpp), qui sont les tops de synthese.
#pragma HLS DATAFLOW

    constexpr int W8 = 8 * PAR_FACTOR;                  // mot ASCII
    constexpr int W2 = 2 * PAR_FACTOR;                  // mot 2 bits
    constexpr int WS = 2 * SMER_SIZE * PAR_FACTOR;      // mot s-mers/hashes/minima

    hls::stream<ap_uint<W8>> s_base_raw;   hls::stream<ap_uint<PAR_FACTOR>> s_base_raw_v;   // m1->m2
    hls::stream<ap_uint<W2>> s_base_2b;    hls::stream<ap_uint<PAR_FACTOR>> s_base_2b_v;    // m2->m3
    hls::stream<ap_uint<W2>> s_base_al;    hls::stream<ap_uint<PAR_FACTOR>> s_base_al_v;    // m3->m4
    hls::stream<ap_uint<WS>> s_smer;       hls::stream<ap_uint<PAR_FACTOR>> s_smer_v;       // m4->m5
    hls::stream<ap_uint<WS>> s_hash;       hls::stream<ap_uint<PAR_FACTOR>> s_hash_v;       // m5->m6
    hls::stream<ap_uint<WS>> s_win;        hls::stream<ap_uint<PAR_FACTOR>> s_win_v;        // m6->m7
    hls::stream<ap_uint<WS>> s_min;        hls::stream<ap_uint<PAR_FACTOR>> s_min_v;        // m7->m8
    hls::stream<ap_uint<WS>> s_ded;        hls::stream<ap_uint<PAR_FACTOR>> s_ded_v;        // m8->m9

#pragma HLS STREAM variable=s_base_raw   depth=2
#pragma HLS STREAM variable=s_base_raw_v depth=2
#pragma HLS STREAM variable=s_base_2b    depth=2
#pragma HLS STREAM variable=s_base_2b_v  depth=2
#pragma HLS STREAM variable=s_base_al    depth=2
#pragma HLS STREAM variable=s_base_al_v  depth=2
#pragma HLS STREAM variable=s_smer       depth=2
#pragma HLS STREAM variable=s_smer_v     depth=2
#pragma HLS STREAM variable=s_hash       depth=2
#pragma HLS STREAM variable=s_hash_v     depth=2
#pragma HLS STREAM variable=s_win        depth=2
#pragma HLS STREAM variable=s_win_v      depth=2
#pragma HLS STREAM variable=s_min        depth=2
#pragma HLS STREAM variable=s_min_v      depth=2
#pragma HLS STREAM variable=s_ded        depth=2
#pragma HLS STREAM variable=s_ded_v      depth=2

    thr_reader<PAR_FACTOR>         (packed_sequence, n_bases, s_base_raw, s_base_raw_v);  // m1
    thr_adapter_hls<PAR_FACTOR>    (s_base_raw, s_base_raw_v, s_base_2b, s_base_2b_v);    // m2
    m3_thr_base_adapter<PAR_FACTOR>(s_base_2b, s_base_2b_v, s_base_al, s_base_al_v);      // m3
    thr_smer_gen<PAR_FACTOR>       (s_base_al, s_base_al_v, s_smer, s_smer_v);            // m4
    thr_hash<PAR_FACTOR>           (s_smer, s_smer_v, s_hash, s_hash_v);                  // m5
    thr_adapter_smer<PAR_FACTOR>   (s_hash, s_hash_v, s_win, s_win_v);                    // m6
    thr_min_v8<PAR_FACTOR>         (s_win, s_win_v, s_min, s_min_v);                      // m7
    thr_dedup_v8<PAR_FACTOR>       (s_min, s_min_v, s_ded, s_ded_v);                      // m8
    thr_store_burst<PAR_FACTOR>    (s_ded, s_ded_v, tab_hash, nMinizrs);                  // m9
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
