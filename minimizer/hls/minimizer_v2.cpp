#include <vector>
#include <cstdint>
#include "ap_int.h"

#define SMER_SIZE         38
#define SMERS_PER_CYCLE   8
#define DEDUP_WINDOW_SIZE 16

// ============================================================
// Hash function
// ============================================================

template <int resu_size>
inline ap_uint<resu_size> hash_u64(uint64_t key) {
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return (ap_uint<resu_size>)key;
}

// ============================================================
// Encode
// ============================================================

static inline uint8_t nucl_encode(uint8_t nucl) {
    return (uint8_t)((nucl >> 1) & 0x03);
}

static inline ap_uint<SMER_SIZE> min_u(
    ap_uint<SMER_SIZE> a,
    ap_uint<SMER_SIZE> b
) {
    return (a < b) ? a : b;
}

// ============================================================
// Reader
// ============================================================

static void reader(
    const ap_uint<64>* packed_sequence,
    uint64_t n_bases,
    std::vector<uint8_t>& out
) {
    out.clear();
    out.reserve(n_bases);

    uint64_t n_words = (n_bases + 7) >> 3;

    for (uint64_t w = 0; w < n_words; w++) {
        uint64_t word = (uint64_t)packed_sequence[w];
        for (int j = 0; j < 8; j++) {
            uint64_t idx = (w << 3) + j;
            if (idx < n_bases)
                out.push_back((word >> (8*j)) & 0xFF);
        }
    }
}

// ============================================================
// Pack (8 bases → 24 bits)
// ============================================================

static void pack(
    const std::vector<uint8_t>& bases,
    std::vector<uint32_t>& out
) {
    out.clear();
    out.reserve((bases.size() + 7) / 8);

    for (size_t i = 0; i < bases.size(); i += 8) {
        uint32_t pack = 0;

        for (int k = 0; k < 8; k++) {
            
            size_t idx = i + k;
            if (idx >= bases.size()) break;

            uint8_t enc = nucl_encode(bases[idx]);
            pack |= ((uint32_t)(enc & 0x3)) << (3*k);
            pack |= ((uint32_t)1) << (3*k + 2);
        }
        out.push_back(pack);
    }
}

// ============================================================
// SMER + Hash
// ============================================================

static void smer(
    const std::vector<uint32_t>& packed24,
    std::vector<uint64_t>& out_hash
) {
    out_hash.clear();

    const int smer_bases = SMER_SIZE / 2;
    ap_uint<SMER_SIZE> cur = 0;
    ap_uint<SMER_SIZE> cur_rc = 0;
    uint64_t count = 0;

    for (uint32_t word : packed24) {
        for (int k = 0; k < 8; k++) {

            if (((word >> (3*k + 2)) & 1) == 0)
                break;

            uint8_t c = (word >> (3*k)) & 0x3;

            cur = (cur << 2) | c;
            cur_rc = (cur_rc >> 2)
                   | ((ap_uint<SMER_SIZE>)(0x2 ^ c) << (SMER_SIZE - 2));

            count++;

            if (count >= (uint64_t)smer_bases) {
                ap_uint<SMER_SIZE> vmin = min_u(cur, cur_rc);
                out_hash.push_back(
                    (uint64_t)hash_u64<SMER_SIZE>((uint64_t)vmin)
                );
            }
        }
    }
}

// ============================================================
// Packet
// ============================================================

struct smer_packet {
    uint64_t v[SMERS_PER_CYCLE];
};

// ============================================================
// Accumulate 
// ============================================================

static void accumulate(
    const std::vector<uint64_t>& in,
    std::vector<smer_packet>& out
) {
    out.clear();
    out.reserve((in.size() + SMERS_PER_CYCLE - 1) / SMERS_PER_CYCLE);

    size_t idx = 0;
    while (idx < in.size()) {

        smer_packet pkt{};
        for (int i = 0; i < SMERS_PER_CYCLE; i++) {
            if (idx < in.size())
                pkt.v[i] = in[idx++];
            else
                pkt.v[i] = 0;   // padding SAFE
        }
        out.push_back(pkt);
    }
}

// ============================================================
// Dedup / Minimizer
// ============================================================

static void dedup(
    const std::vector<smer_packet>& packets,
    size_t n_smers,
    std::vector<uint64_t>& out
) {
    out.clear();

    if (n_smers < DEDUP_WINDOW_SIZE)
        return;

    std::vector<uint64_t> flat;
    flat.reserve(n_smers);

    for (const auto& p : packets) {
        for (int i = 0; i < SMERS_PER_CYCLE; i++) {
            if (flat.size() < n_smers)
                flat.push_back(p.v[i]);
        }
    }

    const size_t limit = flat.size() - DEDUP_WINDOW_SIZE + 1;

    for (size_t i = 0; i < limit; i++) {
        uint64_t m = flat[i];
        for (int j = 1; j < DEDUP_WINDOW_SIZE; j++) {
            if (flat[i + j] < m)
                m = flat[i + j];
        }

        if (out.empty() || out.back() != m)
            out.push_back(m);
    }
}

// ============================================================
// Store
// ============================================================

static void store(
    const std::vector<uint64_t>& in,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMin
) {
    for (size_t i = 0; i < in.size(); i++)
        tab_hash[i] = in[i];

    *nMin = (ap_uint<64>)in.size();
}

// ============================================================
// TOP
// ============================================================

void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n_bases,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
    std::vector<uint8_t> bases;
    std::vector<uint32_t> packed24;
    std::vector<uint64_t> smers;
    std::vector<smer_packet> packets;
    std::vector<uint64_t> minimizers;

    uint64_t n_smers =
        (n_bases >= (SMER_SIZE / 2)) ?
        (uint64_t)n_bases - (SMER_SIZE / 2) + 1 :
        0;

    reader(packed_sequence, (uint64_t)n_bases, bases);
    pack(bases, packed24);
    smer(packed24, smers);
    accumulate(smers, packets);
    dedup(packets, n_smers, minimizers);
    store(minimizers, tab_hash, nMinizrs);
}
