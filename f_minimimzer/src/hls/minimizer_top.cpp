#include "./minimizer_top.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  Tops de synthese concrets. Chaque wrapper porte les interfaces AXI et appelle
//  l'instanciation minimizer<N> correspondante. Choisir le wrapper a synthetiser
//  selon le facteur de parallelisme voulu.
//
//
#define MINIMIZER_WRAPPER(NAME, PF)                                                       \
void NAME(                                                                                \
    const ap_uint< 8 * PF >* packed_sequence,                                             \
    ap_uint<64>*             tab_hash,                                                     \
    ap_uint<64>*             nMinizrs,                                                     \
    ap_uint<64>              n_bases                                                       \
){                                                                                        \
    _Pragma("HLS INTERFACE m_axi     port=packed_sequence offset=slave bundle=gmem0")     \
    _Pragma("HLS INTERFACE m_axi     port=tab_hash        offset=slave bundle=gmem1")     \
    _Pragma("HLS INTERFACE s_axilite port=nMinizrs")                                      \
    _Pragma("HLS INTERFACE s_axilite port=n_bases")                                       \
    _Pragma("HLS INTERFACE s_axilite port=return")                                        \
    minimizer<PF>(packed_sequence, tab_hash, nMinizrs, n_bases);                          \
}
//
MINIMIZER_WRAPPER(minimizer_pf8,   8)
MINIMIZER_WRAPPER(minimizer_pf16, 16)
MINIMIZER_WRAPPER(minimizer_pf32, 32)
MINIMIZER_WRAPPER(minimizer_pf64, 64)
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
