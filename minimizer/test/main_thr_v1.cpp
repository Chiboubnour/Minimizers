#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring> 
#include <fstream> 
#include <inttypes.h>

#include "ap_int.h"

extern int t_thread_minimizer(
    const uint64_t* packed_sequence,
    const int s,
    const int w,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
);

int main(int argc, char * argv[])
{
    printf("\n");
    printf("=============================================\n");
    printf(" Minimizer test   : %s\n", __FILE__);
    printf(" Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("=============================================\n\n");

    std::string i_file;                 
    std::string o_file = "result_thr_v1.txt";

    int s = 19;
    int w = 16;
    bool verbose = true;

    // -----------------------------
    // Parsing des arguments
    // -----------------------------
    for (int i = 1; i < argc; i++)
    {
        std::string argvi = argv[i];

        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i") {
            i_file = argv[++i];
        }
        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o") {
            o_file = argv[++i];
        }
        else if (argvi == "--window" || argvi == "-w") {
            w = std::atoi(argv[++i]);
        }
        else if (argvi == "--smer" || argvi == "-s") {
            s = std::atoi(argv[++i]);
        }
        else if (argvi == "--no-verbose") {
            verbose = false;
        }
        else if (argvi == "-verbose") {
            verbose = true;
        }
        else {
            printf("(EE) Unknown argument: %s\n", argv[i]);
            return -1;
        }
    }

    if (i_file.empty()) {
        printf("(EE) No input file provided. Use -i <file>\n");
        return -1;
    }

    // -----------------------------
    // Lecture de la séquence
    // -----------------------------
    std::string nucl;
    std::ifstream ii(i_file);

    if (ii.is_open()) {
        std::getline(ii, nucl);
        ii.close();
    } else {
        printf("(EE) Impossible to open input file (%s)\n", i_file.c_str());
        return -1;
    }

    const int n = nucl.size();

    if (n == 0) {
        printf("(EE) No data were loaded !\n");
        return -1;
    }

    const int n_smers = n - s + 1;
    const int n_minzr_est = (2 * n_smers) / (w + 1);

    printf("------------ Parameters --------------------\n");
    printf(" Sequence length      : %d\n", n);
    printf(" s-mer size (s)       : %d\n", s);
    printf(" Window size (w)      : %d\n", w);
    printf(" Number of s-mers     : %d\n", n_smers);
    printf(" Estimated minimizers : %d\n", n_minzr_est);
    printf("--------------------------------------------\n\n");

    // -------------------------
    // Pack ASCII sequence en uint64_t
    // -------------------------
    int n64 = (n + 7) / 8;
    std::vector<uint64_t> packed_seq(n64, 0);

    for (int i = 0; i < n; i++) {
        int word_idx = i / 8;
        int byte_idx = i % 8;
        packed_seq[word_idx] |= ((uint64_t)nucl[i] & 0xFF) << (8 * byte_idx);
    }

    ap_uint<64>* tab_hash = new ap_uint<64>[n];
    ap_uint<64> nMin_hw = 0;

    // -------------------------
    // Appel accélérateur
    // -------------------------
    const int nHashs = t_thread_minimizer(
        packed_seq.data(),
        s,
        w,
        tab_hash,
        &nMin_hw
    );

    // -------------------------
    // Affichage résultat
    // -------------------------
    int nPrint = 5;

    printf("\n✅ REAL number of minimizers found : %d\n", nHashs);

    printf(">>> 5 FIRST HASH VALUES\n");
    for (int i = 0; i < std::min(nHashs, nPrint); ++i) {
        printf("  [%4d]  0x%016llX\n", i, tab_hash[i].to_uint64());
    }

    printf("\n>>> 5 LAST HASH VALUES\n");
    int start = (nHashs > nPrint) ? nHashs - nPrint : 0;
    for (int i = start; i < nHashs; ++i) {
        printf("  [%4d]  0x%016llX\n", i, tab_hash[i].to_uint64());
    }

    printf("\n");


    FILE* ff = fopen(o_file.c_str(), "w");
    if (ff) {
        for (int i = 0; i < nHashs; ++i) {
            fprintf(ff, "[%4d] 0x%016llX\n", i, tab_hash[i].to_uint64());
        }
        fclose(ff);
        printf(" Results written in : %s\n", o_file.c_str());
    } else {
        printf("(EE) Impossible to open output file (%s)\n", o_file.c_str());
    }

    delete[] tab_hash;

    printf("\n✅ DONE SUCCESSFULLY\n");

    return 0;
}
