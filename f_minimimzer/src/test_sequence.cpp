//
//  Tests de reference (golden vectors) : des sequences ADN concretes et leurs
//  listes de minimizers attendues (issues de l'implementation de reference).
//  On execute le pipeline complet minimizer<PF>() et on compare la sortie.
//
//  On verifie aussi que le resultat est IDENTIQUE pour tous les PAR_FACTOR
//  (8, 16, 32, 64) : le facteur de parallelisme ne change que le debit.
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
/////////////////////////////////////////////////////////////////////////////////
//
//  Jeux de reference.
//
struct golden_t {
    const char*           seq;
    std::vector<uint64_t> minimizers;
};

static const golden_t GOLDENS[] = {
    {
        "AGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT",
        { 0x00000000719955ddull, 0x00000001cda65704ull, 0x00000004c81719d6ull,
          0x00000002ee272069ull, 0x0000000ea4eb1672ull, 0x0000000c1885d273ull,
          0x000000076423f16eull }
    },
    {
        "GTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGAC",
        { 0x000000014eba320full, 0x000000042ea67e34ull, 0x00000001acb5b83cull,
          0x000000037aa06f13ull, 0x00000003efbf59c7ull, 0x000000018cb17ce9ull,
          0x00000002cb7a4d67ull }
    },
    {
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTG"
        "GTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGAC",
        { 0x0000000cb28c4614ull, 0x0000000135559efeull, 0x000000071172c330ull,
          0x00000009a652c9ffull, 0x00000009e2961f07ull, 0x0000000aa4e82cdfull,
          0x000000035678814full, 0x0000000accc73d13ull, 0x0000000188e7c2c2ull,
          0x000000014eba320full, 0x000000042ea67e34ull, 0x00000001acb5b83cull,
          0x000000037aa06f13ull, 0x00000003efbf59c7ull, 0x000000018cb17ce9ull,
          0x00000002cb7a4d67ull }
    },
};
static constexpr int N_GOLDENS = sizeof(GOLDENS) / sizeof(GOLDENS[0]);
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  Execute le pipeline complet pour un PAR_FACTOR donne et renvoie les minimizers.
//
template<int PF>
static std::vector<uint64_t> run(const std::string& seq)
{
    const uint64_t n_bases = seq.size();
    const uint64_t n_words = (n_bases + PF - 1) / PF;

    std::vector<ap_uint<8 * PF>> mem(n_words, 0);
    for (uint64_t j = 0; j < n_bases; j++)
        mem[j / PF].range(8 * (j % PF) + 7, 8 * (j % PF)) = (uint8_t)seq[j];

    std::vector<ap_uint<64>> tab(1024, 0);
    ap_uint<64>              nMinizrs = 0;

    minimizer<PF>(mem.data(), tab.data(), &nMinizrs, n_bases);

    std::vector<uint64_t> out;
    for (uint64_t i = 0; i < (uint64_t)nMinizrs; i++) out.push_back((uint64_t)tab[i]);
    return out;
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
TEST_CASE("sequence golden : minimizers attendus (PF=8)", "[seq][golden]")
{
    for (int g = 0; g < N_GOLDENS; g++) {
        INFO("jeu golden #" << g);
        const std::vector<uint64_t> got = run<8>(GOLDENS[g].seq);
        const std::vector<uint64_t>& exp = GOLDENS[g].minimizers;

        REQUIRE( got.size() == exp.size() );
        for (size_t i = 0; i < exp.size(); i++) {
            INFO("Minimizer [" << i << "]");
            REQUIRE( got[i] == exp[i] );
        }
    }
}

TEST_CASE("sequence golden : resultat independant de PAR_FACTOR", "[seq][golden]")
{
    for (int g = 0; g < N_GOLDENS; g++) {
        INFO("jeu golden #" << g);
        const std::string         seq = GOLDENS[g].seq;
        const std::vector<uint64_t> g8  = run<8> (seq);
        const std::vector<uint64_t> g16 = run<16>(seq);
        const std::vector<uint64_t> g32 = run<32>(seq);
        const std::vector<uint64_t> g64 = run<64>(seq);
        REQUIRE( g16 == g8 );
        REQUIRE( g32 == g8 );
        REQUIRE( g64 == g8 );
    }
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
