#include <cstdint>
#include <vector>
#include <algorithm>
#include "ref_minimizer.cpp"
#include <omp.h>


template<int S, int W>
static std::vector<uint64_t> process_chunk(
    const char* seq,
    int start,
    int len,
    int n_total,
    int skip_windows
) {
    constexpr int SMER_BITS = 2 * S;
    constexpr uint64_t SMER_MASK =
        (SMER_BITS == 64) ? ~0ULL : ((1ULL << SMER_BITS) - 1);

    int end = std::min(start + len, n_total);
    int chunk_len = end - start;

    if (chunk_len < S)
        return {};

    std::vector<uint64_t> vhash;
    vhash.reserve(chunk_len);

    uint64_t fwd = 0;
    uint64_t rev = 0;

    for (int i = 0; i < chunk_len; ++i) {
        uint64_t b = nucl_encode(seq[start + i]);

        fwd = ((fwd << 2) | b) & SMER_MASK;
        rev = (rev >> 2) | ((0x2ULL ^ b) << (SMER_BITS - 2));

        if (i + 1 >= S) {
            uint64_t canon = std::min(fwd, rev);
            vhash.push_back(hash_u64<S>(canon) & SMER_MASK);
        }
    }

    std::vector<uint64_t> result;

    if ((int)vhash.size() < W)
        return result;

    result.reserve(vhash.size());

    for (size_t i = skip_windows; i + W <= vhash.size(); ++i) {

        uint64_t minv = vhash[i];

        for (int j = 1; j < W; ++j)
            minv = std::min(minv, vhash[i + j]);

        result.push_back(minv);
    }

    return result;
}

template<int S, int W>
int t_minimizer_sequence_omp(
    const char* seq,
    int         n_bases,
    uint64_t*   out_hashes,
    int         n_threads = 0
) {
    constexpr int OVERLAP_BASES = S + W - 2;    // 33
    constexpr int SKIP          = W - 1;       // 15

    if (n_threads <= 0)
        n_threads = omp_get_max_threads();
        
    if (n_threads == 1)
        return t_minimizer_sequence<S,W>(seq, n_bases, out_hashes);
    int chunk_base = std::max(
        (n_bases + n_threads - 1) / n_threads,
        S + W
    );

    std::vector<std::vector<uint64_t>> partial(n_threads);

    #pragma omp parallel for schedule(static) num_threads(n_threads)
    for (int t = 0; t < n_threads; ++t) {

        int start = t * chunk_base;
        if (start >= n_bases) continue;

        int start_ov = (t == 0) ? 0 : (start - OVERLAP_BASES);
        int start_next = (t + 1) * chunk_base;
        int end_ov = (t == n_threads - 1)
                     ? n_bases
                     : std::min(start_next + OVERLAP_BASES, n_bases);

        int len_ov = end_ov - start_ov;
        int skip   = (t == 0) ? 0 : SKIP;

        partial[t] = process_chunk<S, W>(
            seq,
            start_ov,
            len_ov,
            n_bases,
            skip
        );

        if (t < n_threads - 1) {
            int max_owned_windows = start_next - start_ov - (S - 1);
            int owned_start       = skip;
            int max_emit          = std::max(0, max_owned_windows - owned_start);
            if ((int)partial[t].size() > max_emit)
                partial[t].resize(max_emit);
        }
    }

    int out_cnt = 0;
    uint64_t last = UINT64_MAX;
    
    for (int t = 0; t < n_threads; ++t) {
    
        for (uint64_t h : partial[t]) {
    
            if (out_cnt == 0 || h != last) {
                out_hashes[out_cnt++] = h;
                last = h;
            }
        }
    }

    return out_cnt;
}

template int t_minimizer_sequence_omp<19, 16>(
    const char*, int, uint64_t*, int
);