//
//  Testbench du module m9 (store burst), template sur PAR_FACTOR.
//
//  Chaine m1<PF> -> ... -> m9<PF>. m9 ecrit en memoire (dense) les minimizers
//  valides. tab_hash[0..nElem-1] == deduplication consecutive des minima.
//  PF = 8, 16, 32, 64.
//
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
//
#include "./catch2v3/catch_amalgamated.hpp"
//
#include "./hls/m1_thr_reader.hpp"
#include "./hls/m2_thr_base_2b_conv.hpp"
#include "./hls/m3_thr_base_adapter.hpp"
#include "./hls/m4_thr_smer_gen.hpp"
#include "./hls/m5_thr_hash.hpp"
#include "./hls/m6_thr_adapter_smer.hpp"
#include "./hls/m7_thr_min_v8.hpp"
#include "./hls/m8_thr_dedup_v8.hpp"
#include "./hls/m9_thr_store_burst.hpp"
//
//
static const std::string GOLDEN =
    "ACTGACTG" "GTCAGTCA" "AACCTTGG" "GGTTCCAA"
    "ACACACAC" "GTGTGTGT" "TTTTAAAA" "CGCGCGCG";
//
constexpr int SMER_BITS = 2 * SMER_SIZE;
//
static uint8_t ref_code(char b){ switch(b){case 'A':return 0;case 'C':return 1;case 'T':return 2;case 'G':return 3;default:return 0xFF;} }
static uint64_t ref_smer(int k)
{
    uint64_t fwd = 0;
    for (int t = 0; t < SMER_SIZE; t++) fwd |= (uint64_t)ref_code(GOLDEN[k - (SMER_SIZE - 1) + t]) << (2 * (SMER_SIZE - 1 - t));
    uint64_t rev = 0;
    for (int t = 0; t < SMER_SIZE; t++) { const uint64_t base = (fwd >> (2 * t)) & 3; rev |= ((base ^ 2) & 3) << (2 * (SMER_SIZE - 1 - t)); }
    return (fwd < rev) ? fwd : rev;
}
static uint64_t ref_hash(uint64_t key)
{
    key = (~key + (key << 21)); key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8)); key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4)); key ^= key >> 28;
    key = (key + (key << 31));
    const uint64_t mask = (SMER_BITS >= 64) ? ~0ull : ((1ull << SMER_BITS) - 1);
    return key & mask;
}
static uint64_t ref_value(int j) { return ref_hash(ref_smer((SMER_SIZE - 1) + j)); }
static uint64_t ref_min(int t)
{
    uint64_t m = ref_value(t - (WINDOW_SIZE - 1));
    for (int w = 1; w < WINDOW_SIZE; w++) { const uint64_t v = ref_value(t - (WINDOW_SIZE - 1) + w); if (v < m) m = v; }
    return m;
}
static std::vector<uint64_t> expected_minimizers(uint64_t n_bases)
{
    const int H = (int)n_bases - (SMER_SIZE - 1);
    std::vector<uint64_t> v;
    for (int t = WINDOW_SIZE - 1; t < H; t++) { const uint64_t m = ref_min(t); if (v.empty() || m != v.back()) v.push_back(m); }
    return v;
}
template<int PF>
static std::vector<ap_uint<8 * PF>> pack_words(uint64_t n_bases)
{
    const uint64_t n_words = (n_bases + PF - 1) / PF;
    std::vector<ap_uint<8 * PF>> mem(n_words, 0);
    for (uint64_t j = 0; j < n_bases; j++) mem[j / PF].range(8 * (j % PF) + 7, 8 * (j % PF)) = (uint8_t)GOLDEN[j];
    return mem;
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
template<int PF>
static void check(uint64_t n_bases)
{
    INFO("PAR_FACTOR = " << PF << ", n_bases = " << n_bases);

    const std::vector<ap_uint<8 * PF>> mem = pack_words<PF>(n_bases);

    hls::stream<ap_uint<8 * PF>>             s_word;  hls::stream<ap_uint<PF>> s_valid;
    hls::stream<ap_uint<2 * PF>>             s_enc;   hls::stream<ap_uint<PF>> s_enc_v;
    hls::stream<ap_uint<2 * PF>>             s_align; hls::stream<ap_uint<PF>> s_align_v;
    hls::stream<ap_uint<2 * SMER_SIZE * PF>> s_smer;  hls::stream<ap_uint<PF>> s_smer_v;
    hls::stream<ap_uint<2 * SMER_SIZE * PF>> s_hash;  hls::stream<ap_uint<PF>> s_hash_v;
    hls::stream<ap_uint<2 * SMER_SIZE * PF>> s_win;   hls::stream<ap_uint<PF>> s_win_v;
    hls::stream<ap_uint<2 * SMER_SIZE * PF>> s_min;   hls::stream<ap_uint<PF>> s_min_v;
    hls::stream<ap_uint<2 * SMER_SIZE * PF>> s_ded;   hls::stream<ap_uint<PF>> s_ded_v;

    std::vector<ap_uint<64>> tab(256, 0);
    ap_uint<64>              nElem = 0;

    thr_reader<PF>         (mem.data(), n_bases, s_word, s_valid);
    thr_adapter_hls<PF>    (s_word, s_valid, s_enc, s_enc_v);
    m3_thr_base_adapter<PF>(s_enc, s_enc_v, s_align, s_align_v);
    thr_smer_gen<PF>       (s_align, s_align_v, s_smer, s_smer_v);
    thr_hash<PF>           (s_smer, s_smer_v, s_hash, s_hash_v);
    thr_adapter_smer<PF>   (s_hash, s_hash_v, s_win, s_win_v);
    thr_min_v8<PF>         (s_win, s_win_v, s_min, s_min_v);
    thr_dedup_v8<PF>       (s_min, s_min_v, s_ded, s_ded_v);
    thr_store_burst<PF>    (s_ded, s_ded_v, tab.data(), &nElem);

    const std::vector<uint64_t> exp = expected_minimizers(n_bases);
    REQUIRE( (uint64_t)nElem == exp.size() );
    for (size_t i = 0; i < exp.size(); i++) { INFO("tab index i = " << i); REQUIRE( (uint64_t)tab[i] == exp[i] ); }
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
TEST_CASE("m9 store PF=8",  "[m9][store]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<8> (n); }
TEST_CASE("m9 store PF=16", "[m9][store]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<16>(n); }
TEST_CASE("m9 store PF=32", "[m9][store]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<32>(n); }
TEST_CASE("m9 store PF=64", "[m9][store]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<64>(n); }
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
