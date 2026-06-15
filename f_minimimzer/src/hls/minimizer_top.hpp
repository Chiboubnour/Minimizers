#pragma once
#include "./minimizer.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  Tops de synthese concrets : un wrapper par facteur de parallelisme.
//
//  Chacun est une fonction non-template a nom fixe (utilisable comme "top" dans
//  Vitis HLS), portant les interfaces AXI, et appelant minimizer<N>. Le mot
//  d'entree fait 8*N bits (N bases ASCII / mot).
//
//      minimizer_pf8   : PAR_FACTOR =  8   (mot d'entree =  64 bits)
//      minimizer_pf16  : PAR_FACTOR = 16   (mot d'entree = 128 bits)
//      minimizer_pf32  : PAR_FACTOR = 32   (mot d'entree = 256 bits)
//      minimizer_pf64  : PAR_FACTOR = 64   (mot d'entree = 512 bits)
//
extern void minimizer_pf8 (const ap_uint< 64>* packed_sequence, ap_uint<64>* tab_hash, ap_uint<64>* nMinizrs, ap_uint<64> n_bases);
extern void minimizer_pf16(const ap_uint<128>* packed_sequence, ap_uint<64>* tab_hash, ap_uint<64>* nMinizrs, ap_uint<64> n_bases);
extern void minimizer_pf32(const ap_uint<256>* packed_sequence, ap_uint<64>* tab_hash, ap_uint<64>* nMinizrs, ap_uint<64> n_bases);
extern void minimizer_pf64(const ap_uint<512>* packed_sequence, ap_uint<64>* tab_hash, ap_uint<64>* nMinizrs, ap_uint<64> n_bases);
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
