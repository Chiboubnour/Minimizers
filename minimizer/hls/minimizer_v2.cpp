// minimizer_v2.cpp  (VERSION SÛRE - SANS SENTINEL + ANTI SEGFAULT)
#include <vector>
#include <cstdint>
#include "ap_int.h"

#define SMER_SIZE         38
#define SMERS_PER_CYCLE   8
#define DEDUP_WINDOW_SIZE 16

template <int resu_size>
inline ap_uint<resu_size> hash_u64_tpl(uint64_t key) {
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return (ap_uint<resu_size>)key;
}

static inline uint8_t nucl_encode_v1_from_byte(uint8_t nucl) {
    return (uint8_t)((nucl >> 1) & 0x03);
}

static inline ap_uint<SMER_SIZE> min_v1_u(ap_uint<SMER_SIZE> a,
                                           ap_uint<SMER_SIZE> b) {
    return (a < b) ? a : b;
}

static void cpu_thread_reader(
    const ap_uint<64>* packed_sequence,
    uint64_t n_bases,
    std::vector<uint8_t>& out_bases
) {
    out_bases.clear();
    out_bases.reserve(n_bases);

    size_t n_words = (n_bases + 7) / 8;

    for (size_t wi = 0; wi < n_words; ++wi) {
        uint64_t word = (uint64_t)packed_sequence[wi];

        for (int j = 0; j < 8; ++j) {
            size_t idx = wi * 8 + j;
            if (idx >= n_bases) return;

            uint8_t c = (uint8_t)((word >> (8 * j)) & 0xFF);
            out_bases.push_back(c);
        }
    }
}

static void cpu_thread_reader_pack(
    const std::vector<uint8_t>& bases,
    std::vector<uint32_t>& out_packed24
) {
    out_packed24.clear();

    for (size_t i = 0; i < bases.size(); i += 8) {
        uint32_t pack = 0;

        for (int j = 0; j < 8; ++j) {
            size_t idx = i + j;
            if (idx >= bases.size()) break;

            uint8_t enc = nucl_encode_v1_from_byte(bases[idx]);

            pack |= ((uint32_t)(enc & 0x3)) << (3 * j);
            pack |= ((uint32_t)1) << (3 * j + 2); // valid bit
        }

        out_packed24.push_back(pack);
    }
}

static void cpu_thread_smer(
    const std::vector<uint32_t>& packed24,
    std::vector<uint64_t>& out_vhash
) {
    out_vhash.clear();

    const int smer_bases = SMER_SIZE / 2;

    ap_uint<SMER_SIZE> current = 0;
    ap_uint<SMER_SIZE> cur_inv = 0;
    uint64_t base_count = 0;

    for (uint32_t word : packed24) {
        for (int k = 0; k < SMERS_PER_CYCLE; ++k) {
            uint32_t trio = (word >> (3 * k)) & 0x7;
            if (((trio >> 2) & 1) == 0) break;

            uint8_t enc = trio & 0x3;

            current = (current << 2) | enc;

            ap_uint<SMER_SIZE> inv_part = (0x2 ^ enc) & 0x3;
            cur_inv = (cur_inv >> 2) | (inv_part << (SMER_SIZE - 2));

            base_count++;

            if (base_count >= smer_bases) {
                ap_uint<SMER_SIZE> vmin = min_v1_u(current, cur_inv);
                ap_uint<SMER_SIZE> h = hash_u64_tpl<SMER_SIZE>((uint64_t)vmin);
                out_vhash.push_back((uint64_t)h);
            }
        }
    }
}

static void cpu_thread_dedup(
    const std::vector<uint64_t>& vhash,
    std::vector<uint64_t>& out_minimizers,
    int window
) {
    out_minimizers.clear();
    if (vhash.size() < (size_t)window) return;

    for (size_t i = 0; i + window <= vhash.size(); ++i) {
        uint64_t minv = vhash[i];

        for (int j = 1; j < window; ++j) {
            if (vhash[i + j] < minv)
                minv = vhash[i + j];
        }

        if (out_minimizers.empty() || out_minimizers.back() != minv) {
            out_minimizers.push_back(minv);
        }
    }
}

static void cpu_thread_store(
    const std::vector<uint64_t>& minimizers,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
    for (size_t i = 0; i < minimizers.size(); ++i) {
        tab_hash[i] = minimizers[i];
    }

    if (nMinizrs)
        *nMinizrs = minimizers.size();
}

 void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n_bases,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
    std::vector<uint8_t> bases;
    std::vector<uint32_t> packed24;
    std::vector<uint64_t> vhash;
    std::vector<uint64_t> minimizers;

    cpu_thread_reader(packed_sequence, (uint64_t)n_bases, bases);
    cpu_thread_reader_pack(bases, packed24);
    cpu_thread_smer(packed24, vhash);
    cpu_thread_dedup(vhash, minimizers, DEDUP_WINDOW_SIZE);
    cpu_thread_store(minimizers, tab_hash, nMinizrs);
}
