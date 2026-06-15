#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m1 : lecteur memoire (template sur PAR_FACTOR).
//
//  Un mot memoire fait 8*PAR_FACTOR bits = PAR_FACTOR bases ASCII. Chaque mot est
//  transmis tel quel avec un masque valid de PAR_FACTOR bits (les r bits de poids
//  faible si seulement r bases restent, sinon tout a 1). Un paquet terminal (0,0)
//  clot les flux.
//
//  NB : module header-only (template) -> definition dans le .hpp.
//
//
//  Masque des "reste" bases valides (sur PAR_FACTOR bits).
//
template<int PAR_FACTOR>
inline ap_uint<PAR_FACTOR> valid_bits(const ap_uint<64> reste)
{
    if (reste >= PAR_FACTOR) return ~ap_uint<PAR_FACTOR>(0);          // PAR_FACTOR bases valides
    return (ap_uint<PAR_FACTOR>(1) << reste) - 1;                     // r bases valides
}
//
//
//
template<int PAR_FACTOR>
void thr_reader(
    const ap_uint< 8 * PAR_FACTOR >*      base_ptr_i,
    const ap_uint<            64   >      n_bases_i,
    hls::stream< ap_uint< 8 * PAR_FACTOR > >& base_stream_o,
    hls::stream< ap_uint<     PAR_FACTOR > >& base_valid_o
){
#pragma HLS INLINE off
    //
    // chaque mot contient PAR_FACTOR bases (8 bits ASCII chacune)
    //
    ap_uint<64> n_words  = (n_bases_i + PAR_FACTOR - 1) / PAR_FACTOR;
    ap_uint<64> base_cnt =  n_bases_i;
    for (ap_uint<64> w = 0; w < n_words; w++)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint< 8 * PAR_FACTOR > word  = base_ptr_i[w];
        const ap_uint<     PAR_FACTOR > valid = valid_bits<PAR_FACTOR>( base_cnt );

        base_stream_o.write(word);
        base_valid_o.write (valid);

        base_cnt -= PAR_FACTOR;
    }
    //
    // paquet terminal
    //
    base_stream_o.write(0);
    base_valid_o.write (0);
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
