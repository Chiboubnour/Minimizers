//
//  Testbench unitaire du module m3 (aligneur de s-mer), template sur PAR_FACTOR.
//
//  Chaine m1<PF> -> m2<PF> -> m3<PF>. m3 recadre le flux de bases 2 bits :
//  warm-up de (SMER_SIZE-1) bases livre en first_rounds paquets pleins + 1 paquet
//  de n_last_round bases, puis le reste en paquets pleins. Convention : seul
//  "valid" fait foi (donnees des slots invalides = don't-care). Valide pour
//  PF = 8, 16, 32, 64 (PF=32/64 -> first_rounds=0).
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
//  Cadencage attendu des masques valid de m3 (data don't-care).
//
template<int PF>
static std::vector<ap_uint<PF>> expected_valids(uint64_t n_bases)
{
    constexpr int first_rounds = (SMER_SIZE - 1) / PF;
    constexpr int n_last_round = (SMER_SIZE - 1) % PF;

    std::vector<ap_uint<PF>> v;
    uint64_t idx = 0;
    int      beat = 0;
    while (idx < n_bases) {
        const int target = (beat <  first_rounds) ? PF
                         : (beat == first_rounds) ? n_last_round
                                                  : PF;
        const int take = (target < (int)(n_bases - idx)) ? target : (int)(n_bases - idx);
        v.push_back(((ap_uint<PF>)1 << take) - 1);
        idx += take;
        beat++;
    }
    v.push_back(0);   // terminal
    return v;
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
template<int PF>
struct beats {
    std::vector<ap_uint<2 * PF>> words;
    std::vector<ap_uint<PF>>     valids;
};

template<int PF>
static beats<PF> run_chain(uint64_t n_bases)
{
    const std::vector<ap_uint<8 * PF>> mem = pack_words<PF>(n_bases);

    hls::stream<ap_uint<8 * PF>> s_word;   hls::stream<ap_uint<PF>> s_valid;
    hls::stream<ap_uint<2 * PF>> s_enc;    hls::stream<ap_uint<PF>> s_enc_v;
    hls::stream<ap_uint<2 * PF>> s_align;  hls::stream<ap_uint<PF>> s_align_v;

    thr_reader<PF>         (mem.data(), n_bases, s_word, s_valid);
    thr_adapter_hls<PF>    (s_word, s_valid, s_enc, s_enc_v);
    m3_thr_base_adapter<PF>(s_enc, s_enc_v, s_align, s_align_v);

    beats<PF> out;
    while (!s_align.empty()) {
        out.words.push_back (s_align.read());
        out.valids.push_back(s_align_v.read());
    }
    return out;
}

template<int PF>
static std::string decode_beats(const beats<PF>& out)
{
    std::string seq;
    for (size_t w = 0; w + 1 < out.words.size(); w++) {
        const ap_uint<PF> valid = out.valids[w];
        for (int i = 0; i < PF; i++)
            if (valid[i]) seq.push_back(ref_base((uint8_t)out.words[w].range(2 * i + 1, 2 * i)));
    }
    return seq;
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
template<int PF>
static void check(uint64_t n_bases)
{
    const beats<PF>                got = run_chain<PF>(n_bases);
    const std::vector<ap_uint<PF>> exp = expected_valids<PF>(n_bases);

    INFO("PAR_FACTOR = " << PF << ", n_bases = " << n_bases);

    REQUIRE( got.valids.size() == exp.size() );          // meme nombre de battements
    for (size_t w = 0; w < exp.size(); w++) {
        INFO("battement w = " << w);
        REQUIRE( got.valids[w] == exp[w] );              // seul valid fait foi
    }

    REQUIRE( decode_beats<PF>(got) == GOLDEN.substr(0, n_bases) );   // aucune base perdue
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
TEST_CASE("m3 aligneur PF=8",  "[m3][align]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<8> (n); }
TEST_CASE("m3 aligneur PF=16", "[m3][align]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<16>(n); }
TEST_CASE("m3 aligneur PF=32", "[m3][align]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<32>(n); }
TEST_CASE("m3 aligneur PF=64", "[m3][align]") { const uint64_t n = GENERATE(64,63,62,61,60,59,58,57); check<64>(n); }
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
