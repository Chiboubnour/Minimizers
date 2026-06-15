//
//  Testbench top-level : appelle directement minimizer<PF>() (le pipeline assemble
//  en dataflow) et verifie le tableau de minimizers ecrit en memoire.
//
//  Test d'integration du pipeline complet via son point d'entree unique, pour
//  PAR_FACTOR = 8, 16, 32, 64.
//
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
//
#include "./catch2v3/catch_amalgamated.hpp"
//
#include "./hls/minimizer.hpp"
//
//
static_assert(true, "");
//
constexpr int SMER_BITS = 2 * SMER_SIZE;
//
static const std::string GOLDEN =
    "ACTGACTG" "GTCAGTCA" "AACCTTGG" "GGTTCCAA"
    "ACACACAC" "GTGTGTGT" "TTTTAAAA" "CGCGCGCG";
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

    std::vector<ap_uint<8 * PF>> seq      = pack_words<PF>(n_bases);
    std::vector<ap_uint<64>>     tab(256, 0);
    ap_uint<64>                  nMinizrs = 0;

    minimizer<PF>(seq.data(), tab.data(), &nMinizrs, n_bases);

    const std::vector<uint64_t> exp = expected_minimizers(n_bases);
    REQUIRE( (uint64_t)nMinizrs == exp.size() );
    for (size_t i = 0; i < exp.size(); i++) { INFO("tab index i = " << i); REQUIRE( (uint64_t)tab[i] == exp[i] ); }
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
TEST_CASE("minimizer top PF=8",  "[top]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<8> (n); }
TEST_CASE("minimizer top PF=16", "[top]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<16>(n); }
TEST_CASE("minimizer top PF=32", "[top]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<32>(n); }
TEST_CASE("minimizer top PF=64", "[top]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<64>(n); }
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
