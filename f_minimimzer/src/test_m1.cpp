//
//  Testbench unitaire du module m1 (thr_reader), template sur PAR_FACTOR.
//
//  m1 est desormais parametre par PAR_FACTOR : un mot memoire fait 8*PAR_FACTOR
//  bits (PAR_FACTOR bases ASCII). On valide pour PAR_FACTOR = 8, 16, 32, 64.
//
//  Pour chaque longueur on packe GOLDEN dans des mots de 8*PF bits, on execute
//  thr_reader<PF>, et on verifie : nombre de mots (ceil(n/PF) + terminal), recopie
//  fidele, masque valid (r bits de poids faible), terminal, et reconstruction.
//
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
//
#include "./catch2v3/catch_amalgamated.hpp"
//
#include "./hls/m1_thr_reader.hpp"
//
//
static const std::string GOLDEN =
    "ACTGACTG" "GTCAGTCA" "AACCTTGG" "GGTTCCAA"
    "ACACACAC" "GTGTGTGT" "TTTTAAAA" "CGCGCGCG";   // 64 bases
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  Packe GOLDEN[0..n_bases-1] dans des mots de 8*PF bits (PF bases ASCII / mot).
//
template<int PF>
static std::vector<ap_uint<8 * PF>> pack_words(uint64_t n_bases)
{
    const uint64_t n_words = (n_bases + PF - 1) / PF;
    std::vector<ap_uint<8 * PF>> mem(n_words, 0);
    for (uint64_t j = 0; j < n_bases; j++) {
        const uint64_t w = j / PF;
        const uint64_t i = j % PF;
        mem[w].range(8 * i + 7, 8 * i) = (uint8_t)GOLDEN[j];
    }
    return mem;
}
//
//  Masque valid attendu (PF bits) pour 'reste' bases restantes (oracle independant).
//
template<int PF>
static ap_uint<PF> ref_valid(uint64_t reste)
{
    if (reste >= (uint64_t)PF) return ~ap_uint<PF>(0);
    return (ap_uint<PF>(1) << reste) - 1;
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

    hls::stream<ap_uint<8 * PF>> s_word;
    hls::stream<ap_uint<PF>>     s_valid;

    thr_reader<PF>(mem.data(), n_bases, s_word, s_valid);

    std::vector<ap_uint<8 * PF>> ow;
    std::vector<ap_uint<PF>>     ov;
    while (!s_word.empty()) { ow.push_back(s_word.read()); ov.push_back(s_valid.read()); }

    INFO("PAR_FACTOR = " << PF << ", n_bases = " << n_bases);

    REQUIRE( ow.size() == n_words + 1 );          // n_words + paquet terminal

    uint64_t    remaining = n_bases;
    std::string decoded;
    for (uint64_t w = 0; w < n_words; w++) {
        INFO("word index = " << w);
        REQUIRE( ow[w] == mem[w] );               // recopie fidele du mot
        REQUIRE( ov[w] == ref_valid<PF>(remaining) );
        for (int i = 0; i < PF; i++) {
            if (ov[w][i]) decoded.push_back((char)(uint8_t)ow[w].range(8 * i + 7, 8 * i));
        }
        remaining -= PF;
    }

    // paquet terminal
    REQUIRE( ow.back() == 0 );
    REQUIRE( ov.back() == 0 );

    // reconstruction
    REQUIRE( decoded == GOLDEN.substr(0, n_bases) );
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
TEST_CASE("m1 thr_reader PF=8",  "[m1][reader]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<8> (n); }
TEST_CASE("m1 thr_reader PF=16", "[m1][reader]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<16>(n); }
TEST_CASE("m1 thr_reader PF=32", "[m1][reader]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<32>(n); }
TEST_CASE("m1 thr_reader PF=64", "[m1][reader]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<64>(n); }
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
