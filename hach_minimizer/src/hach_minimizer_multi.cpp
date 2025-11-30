#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <omp.h>


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
extern "C"  void minimizer_cpu_multi(const uint64_t *input_minimizers, uint64_t *output_hashes, unsigned int data_size_words) {
    #pragma omp parallel for num_threads(4)
    for (unsigned int idx = 0; idx < data_size_words; idx++) {
        output_hashes[idx] = bfc_hash_64(input_minimizers[idx], ~0ULL);
    }
}
