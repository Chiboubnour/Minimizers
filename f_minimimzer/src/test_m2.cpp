//
//  Testbench unitaire du module m2 (encodeur 2 bits), template sur PAR_FACTOR.
//
//  Chaine m1<PF> -> m2<PF>. m2 traduit chaque base ASCII en code 2 bits
//  (A=00, C=01, T=10, G=11) et relaie le masque valid. Valide pour PF = 8,16,32,64.
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
//
//
static const std::string GOLDEN =
    "ACTGACTG" "GTCAGTCA" "AACCTTGG" "GGTTCCAA"
    "ACACACAC" "GTGTGTGT" "TTTTAAAA" "CGCGCGCG";   // 64 bases
//
static uint8_t ref_code(char base)
{
    switch (base) {
        case 'A': return 0;  case 'C': return 1;
        case 'T': return 2;  case 'G': return 3;
        default:  return 0xFF;
    }
}
static char ref_base(uint8_t code)
{
    static const char lut[4] = { 'A', 'C', 'T', 'G' };
    return lut[code & 3];
}
template<int PF>
static ap_uint<PF> ref_valid(uint64_t reste)
{
    if (reste >= (uint64_t)PF) return ~ap_uint<PF>(0);
    return (ap_uint<PF>(1) << reste) - 1;
}
//
//  Packe GOLDEN[0..n-1] dans des mots de 8*PF bits (PF bases ASCII / mot).
//
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
template<int PF>
static void check(uint64_t n_bases)
{
    const std::vector<ap_uint<8 * PF>> mem     = pack_words<PF>(n_bases);
    const uint64_t                     n_words = (n_bases + PF - 1) / PF;

    hls::stream<ap_uint<8 * PF>> s_word;  hls::stream<ap_uint<PF>> s_valid;   // m1->m2
    hls::stream<ap_uint<2 * PF>> s_enc;   hls::stream<ap_uint<PF>> s_enc_v;   // m2->out

    thr_reader<PF>     (mem.data(), n_bases, s_word, s_valid);   // m1
    thr_adapter_hls<PF>(s_word, s_valid, s_enc, s_enc_v);        // m2

    std::vector<ap_uint<2 * PF>> ow;  std::vector<ap_uint<PF>> ov;
    while (!s_enc.empty()) { ow.push_back(s_enc.read()); ov.push_back(s_enc_v.read()); }

    INFO("PAR_FACTOR = " << PF << ", n_bases = " << n_bases);

    REQUIRE( ow.size() == n_words + 1 );          // m2 est 1:1 + terminal

    uint64_t    remaining = n_bases;
    std::string decoded;
    for (uint64_t w = 0; w < n_words; w++) {
        INFO("word index = " << w);
        REQUIRE( ov[w] == ref_valid<PF>(remaining) );      // valid relaye inchange
        for (int i = 0; i < PF; i++) {
            if (ov[w][i]) {
                const char    base = GOLDEN[w * PF + i];
                const uint8_t got  = (uint8_t)ow[w].range(2 * i + 1, 2 * i);
                INFO("base i = " << i << " ('" << base << "')");
                REQUIRE( got == ref_code(base) );
                decoded.push_back(ref_base(got));
            }
        }
        remaining -= PF;
    }

    // paquet terminal + reconstruction
    REQUIRE( ow.back() == 0 );
    REQUIRE( ov.back() == 0 );
    REQUIRE( decoded == GOLDEN.substr(0, n_bases) );
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
TEST_CASE("m2 encodeur 2b PF=8",  "[m2][conv]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<8> (n); }
TEST_CASE("m2 encodeur 2b PF=16", "[m2][conv]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<16>(n); }
TEST_CASE("m2 encodeur 2b PF=32", "[m2][conv]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<32>(n); }
TEST_CASE("m2 encodeur 2b PF=64", "[m2][conv]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<64>(n); }
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
