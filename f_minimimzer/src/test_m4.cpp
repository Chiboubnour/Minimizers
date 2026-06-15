//
//  Testbench unitaire du module m4 (s-mers canoniques), template sur PAR_FACTOR.
//
//  Chaine m1<PF> -> ... -> m4<PF>. Verification par sequence aplatie : les s-mers
//  valides en sortie, lus dans l'ordre des positions, doivent valoir
//  ref_smer(k) pour k = SMER_SIZE-1 .. n_bases-1. Valide pour PF = 8,16,32,64.
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
//
//
static const std::string GOLDEN =
    "ACTGACTG" "GTCAGTCA" "AACCTTGG" "GGTTCCAA"
    "ACACACAC" "GTGTGTGT" "TTTTAAAA" "CGCGCGCG";   // 64 bases
//
constexpr int SMER_BITS = 2 * SMER_SIZE;
//
static uint8_t ref_code(char base)
{
    switch (base) {
        case 'A': return 0;  case 'C': return 1;
        case 'T': return 2;  case 'G': return 3;
        default:  return 0xFF;
    }
}
//  s-mer canonique se terminant a la base k
static uint64_t ref_smer(int k)
{
    uint64_t fwd = 0;
    for (int t = 0; t < SMER_SIZE; t++)
        fwd |= (uint64_t)ref_code(GOLDEN[k - (SMER_SIZE - 1) + t]) << (2 * (SMER_SIZE - 1 - t));
    uint64_t rev = 0;
    for (int t = 0; t < SMER_SIZE; t++) {
        const uint64_t base = (fwd >> (2 * t)) & 3;
        rev |= ((base ^ 2) & 3) << (2 * (SMER_SIZE - 1 - t));
    }
    return (fwd < rev) ? fwd : rev;
}
template<int PF>
static std::vector<ap_uint<8 * PF>> pack_words(uint64_t n_bases)
{
    const uint64_t n_words = (n_bases + PF - 1) / PF;
    std::vector<ap_uint<8 * PF>> mem(n_words, 0);
    for (uint64_t j = 0; j < n_bases; j++)
        mem[j / PF].range(8 * (j % PF) + 7, 8 * (j % PF)) = (uint8_t)GOLDEN[j];
    return mem;
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  Execute m1<PF> -> ... -> m4<PF>, aplatit les s-mers valides dans l'ordre.
//
template<int PF>
static std::vector<uint64_t> run_flat(uint64_t n_bases)
{
    const std::vector<ap_uint<8 * PF>> mem = pack_words<PF>(n_bases);

    hls::stream<ap_uint<8 * PF>>            s_word;   hls::stream<ap_uint<PF>> s_valid;
    hls::stream<ap_uint<2 * PF>>            s_enc;    hls::stream<ap_uint<PF>> s_enc_v;
    hls::stream<ap_uint<2 * PF>>            s_align;  hls::stream<ap_uint<PF>> s_align_v;
    hls::stream<ap_uint<2 * SMER_SIZE * PF>> s_smer;  hls::stream<ap_uint<PF>> s_smer_v;

    thr_reader<PF>         (mem.data(), n_bases, s_word, s_valid);
    thr_adapter_hls<PF>    (s_word, s_valid, s_enc, s_enc_v);
    m3_thr_base_adapter<PF>(s_enc, s_enc_v, s_align, s_align_v);
    thr_smer_gen<PF>       (s_align, s_align_v, s_smer, s_smer_v);

    std::vector<uint64_t> flat;
    bool terminal_seen = false;
    while (!s_smer.empty()) {
        const ap_uint<2 * SMER_SIZE * PF> w = s_smer.read();
        const ap_uint<PF>                 v = s_smer_v.read();
        if (v == 0) { terminal_seen = true; continue; }
        for (int i = 0; i < PF; i++)
            if (v[i]) flat.push_back((uint64_t)w.range((i + 1) * SMER_BITS - 1, i * SMER_BITS));
    }
    REQUIRE(terminal_seen);
    return flat;
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
template<int PF>
static void check(uint64_t n_bases)
{
    INFO("PAR_FACTOR = " << PF << ", n_bases = " << n_bases);

    std::vector<uint64_t> exp;
    for (int k = SMER_SIZE - 1; k < (int)n_bases; k++) exp.push_back(ref_smer(k));

    const std::vector<uint64_t> got = run_flat<PF>(n_bases);

    REQUIRE( got.size() == exp.size() );
    for (size_t j = 0; j < exp.size(); j++) {
        INFO("s-mer index j = " << j);
        REQUIRE( got[j] == exp[j] );
    }
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
TEST_CASE("m4 s-mers PF=8",  "[m4][smer]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<8> (n); }
TEST_CASE("m4 s-mers PF=16", "[m4][smer]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<16>(n); }
TEST_CASE("m4 s-mers PF=32", "[m4][smer]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<32>(n); }
TEST_CASE("m4 s-mers PF=64", "[m4][smer]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<64>(n); }
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
