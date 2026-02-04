#include <cstdint>
#include <ap_int.h>
#include <hls_stream.h>

#define HASH_WIDTH 64

inline ap_uint<64> bfc_hash_64(ap_uint<64> key) {
#pragma HLS INLINE
    ap_uint<64> mask = ~0ULL;

    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;

    return key;
}

void minimizer(
    const ap_uint<64>* input_minimizers,
    ap_uint<64>*       output_hashes,
    unsigned int       data_size_words
) {

#ifdef _PRAGMA_HLS_
#pragma HLS INTERFACE m_axi port=input_minimizers offset=slave bundle=gmem_in  max_read_burst_length=256
#pragma HLS INTERFACE m_axi port=output_hashes   offset=slave bundle=gmem_out  max_write_burst_length=256
#pragma HLS INTERFACE s_axilite port=data_size_words
#pragma HLS INTERFACE s_axilite port=return
#endif

    main_loop:
    for (unsigned int i = 0; i < data_size_words; i++) {

#ifdef _PRAGMA_HLS_
        #pragma HLS PIPELINE II=1
#endif

        ap_uint<64> in  = input_minimizers[i];
        ap_uint<64> out = bfc_hash_64(in);
        output_hashes[i] = out;
    }
}
