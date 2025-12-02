#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <inttypes.h>

#include "ap_int.h"

// Attention : le top HLS que je t'ai fourni exporte un wrapper C nommé minimizer_hls_top
extern "C" void minimizer_hls_top(
    const ap_uint<64>* packed_sequence,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
);

int main(int argc, char * argv[])
{
    printf("# \n");
    printf("# \n");
    printf("# Minimizer test   : %s\n", __FILE__);
    printf("# Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("# \n");

    std::string i_file = "./dataset/b64/test_64.nuc";
    std::string o_file = "result.txt";

    int s = 19;
    int w = 16;
    bool verbose = true;

    // parse args (simple)
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--input" || a == "-i") {
            if (i+1 < argc) i_file = argv[++i];
        } else if (a == "--output" || a == "-o") {
            if (i+1 < argc) o_file = argv[++i];
        } else if (a == "--no-verbose") verbose = false;
        else if (a == "-verbose") verbose = true;
        else { printf("(EE) Unknown argument: %s\n", argv[i]); return EXIT_FAILURE; }
    }

    // read first line of the input file (sequence)
    std::string nucl;
    std::ifstream ii(i_file);
    if (!ii.is_open()) {
        printf("(EE) Impossible to open input file (%s)\n", i_file.c_str());
        return EXIT_FAILURE;
    }
    std::getline(ii, nucl);
    ii.close();

    const int n = static_cast<int>(nucl.size());
    if (n <= 0) {
        printf("(EE) No data were loaded !\n");
        return EXIT_FAILURE;
    }

    const int n_smers = n - s + 1;
    const int n_minzr = 2 * n_smers / (w + 1);

    printf("Taille de la sequence     : %d\n", n);
    printf("Nombre de s=%d-mers       : %d\n", s, n_smers);
    printf("Taille de la window.      : %d\n", w);
    printf("Nombre de (%d) minimizers : %d\n", s, n_minzr);

    // ---------------------------------------------------------------------
    // IMPORTANT : le reader HLS lit au maximum 64 mots (64 * 8 = 512 chars) - comme la référence.
    // On construit un tableau de 64 ap_uint<64> et on met les mots restants à 0.
    // ---------------------------------------------------------------------
    const int MAX_WORDS = 64;
    const int n_words = (n + 7) / 8;
    ap_uint<64>* packed_sequence = new ap_uint<64>[MAX_WORDS];
    // init à 0
    for (int i = 0; i < MAX_WORDS; ++i) packed_sequence[i] = 0;

    // pack LSB-first (8 bytes par mot)
    for (int wi = 0; wi < n_words && wi < MAX_WORDS; ++wi) {
        ap_uint<64> w64 = 0;
        for (int j = 0; j < 8; ++j) {
            int idx = wi * 8 + j;
            if (idx < n) {
                uint64_t c = static_cast<uint8_t>(nucl[idx]);
                w64 |= ( (uint64_t)c << (8 * j) );
            }
        }
        packed_sequence[wi] = w64;
    }
    // Si la séquence dépasse 64*8 = 512 caractères, on tronque (comme la référence).
    if (n_words > MAX_WORDS) {
        printf("(WW) Sequence truncated to %d chars for HLS reader (max %d words).\n", MAX_WORDS*8, MAX_WORDS);
    }

    // allocate output
    const int max_out = n; // safe upper bound
    ap_uint<64>* tab_hash = new ap_uint<64>[ (max_out>0) ? max_out : 1 ];
    for (int i = 0; i < max_out; ++i) tab_hash[i] = 0;
    ap_uint<64> nMinizrs = 0;

    // Call HLS top-level (wrapper)
    minimizer_hls_top(packed_sequence, tab_hash, &nMinizrs);

    // convert count
    uint64_t found = nMinizrs.to_uint64();
    printf("HLS returned %llu minimizers\n", (unsigned long long)found);

    // Print results
    if (verbose) {
        for (uint64_t i = 0; i < found; ++i) {
            printf("s-mer [%2llu] : hash 0x%16.16llX\n",
                   (unsigned long long)i, (unsigned long long) tab_hash[i].to_uint64());
        }
    }

    // write to file
    FILE* ff = fopen(o_file.c_str(), "w");
    if (ff) {
        for (uint64_t i = 0; i < found; ++i) {
            fprintf(ff, "s-mer [%2llu] : hash 0x%16.16llX\n",
                    (unsigned long long)i, (unsigned long long) tab_hash[i].to_uint64());
        }
        fclose(ff);
    } else {
        printf("(EE) Impossible d'ouvrir le fichier de sortie (%s)\n", o_file.c_str());
    }

    delete[] tab_hash;
    delete[] packed_sequence;
    return 0;
}
