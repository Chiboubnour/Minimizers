#include <cstdint>

#define HASH_WIDTH 64

inline uint64_t bfc_hash_64(uint64_t key, uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

void parallel_hash_calc(uint64_t input_data, uint64_t& output_data) {
    uint64_t mask = ~0ULL;
    output_data = bfc_hash_64(input_data, mask);
}

void minimizer(
    const uint64_t* input_minimizers,
    uint64_t*       output_hashes,
    unsigned int    data_size_words
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INTERFACE m_axi port=input_minimizers offset=slave bundle=gmem_in max_read_burst_length=256
    #pragma HLS INTERFACE m_axi port=output_hashes offset=slave bundle=gmem_out max_write_burst_length=256
    #pragma HLS INTERFACE s_axilite port=data_size_words
    #pragma HLS INTERFACE s_axilite port=return
#endif

    const unsigned int BURST_SIZE = 1024;

    for (unsigned int burst_idx = 0; burst_idx < data_size_words; burst_idx += BURST_SIZE) {
        unsigned int remaining = data_size_words - burst_idx;
        unsigned int current_burst = (remaining < BURST_SIZE) ? remaining : BURST_SIZE;

        for (unsigned int word_idx = 0; word_idx < current_burst; word_idx++) {
#ifdef _PRAGMA_HLS_
            #pragma HLS PIPELINE II=1
#endif
            unsigned int idx = burst_idx + word_idx;
            uint64_t input_data = input_minimizers[idx];
            uint64_t output_data;

            parallel_hash_calc(input_data, output_data);
            output_hashes[idx] = output_data;
        }
    }
}
