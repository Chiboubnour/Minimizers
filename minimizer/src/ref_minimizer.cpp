#include <cstdint>
#include <vector>
#include <algorithm>


static inline uint64_t nucl_encode(char c) {
    return (c >> 1) & 0x3ULL;
}


template<int S>
inline uint64_t hash_u64(uint64_t key) {
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return key;
}


template<int S, int W>
int t_minimizer_sequence(
    const char* seq,
    int n_bases,
    uint64_t* out_hashes
) {
    constexpr int SMER_BITS = 2 * S;
    constexpr uint64_t SMER_MASK = (SMER_BITS == 64) ? ~0ULL : ((1ULL << SMER_BITS) - 1);

    uint64_t fwd = 0;
    uint64_t rev = 0;

    std::vector<uint64_t> vhash;
    vhash.reserve(n_bases);

    // --- Generate canonical s-mer hashes ---
    for (int i = 0; i < n_bases; ++i) {
        uint64_t b = nucl_encode(seq[i]);

        fwd = ((fwd << 2) | b) & SMER_MASK;
        rev = (rev >> 2) | ((0x2 ^ b) << (SMER_BITS - 2));

        if (i + 1 >= S) {
            uint64_t canon = std::min(fwd, rev);
            uint64_t h = hash_u64<S>(canon) & SMER_MASK; 
            vhash.push_back(h);
        }
    }

    // --- Apply minimizer window ---
    int out_cnt = 0;
    uint64_t last_min = UINT64_MAX;

    for (size_t i = 0; i + W <= vhash.size(); ++i) {
        uint64_t minv = vhash[i];
        for (int j = 1; j < W; ++j)
            minv = std::min(minv, vhash[i + j]);

        if (out_cnt == 0 || minv != last_min) {
            out_hashes[out_cnt++] = minv;
            last_min = minv;
        }
    }

    return out_cnt;
}

template int t_minimizer_sequence<19,16>(
    const char*, int, uint64_t*
);
