#include <ap_int.h>
#include <hls_stream.h>

void reader_hls_64(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n_bases,
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid
) {

    ap_uint<64> n_words = (n_bases + 7) >> 3;

    READ_LOOP:
    for (ap_uint<64> w = 0; w < n_words; w++) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> word = packed_sequence[w];

        ap_uint<8> valid;

        if (w == n_words - 1 && (n_bases & 7)) {
            valid = n_bases & 7;
        } else {
            valid = 8;
        }

        base_stream.write(word);
        base_valid.write(valid);
    }

    base_stream.write(0);
    base_valid.write(0);
}
