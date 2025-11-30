#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>

#define HASH_WIDTH 64
#define DATA_SIZE_BYTES (32 * 1024 * 1024)  
#define NUM_WORDS (DATA_SIZE_BYTES / sizeof(uint64_t))

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


inline void parallel_hash_calc(uint64_t input_data, uint64_t &output_data) {
    uint64_t mask = ~0ULL;
    output_data = bfc_hash_64(input_data, mask);
}

extern "C" void minimizer_cpu_mono(const uint64_t *input_minimizers, uint64_t *output_hashes, unsigned int data_size_words) {
    for (unsigned int idx = 0; idx < data_size_words; idx++) {
        uint64_t input_data = input_minimizers[idx];
        uint64_t output_data;
        parallel_hash_calc(input_data, output_data);
        output_hashes[idx] = output_data;
    }
}

