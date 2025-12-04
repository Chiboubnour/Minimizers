#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstring> 
#include <fstream> 
#include <inttypes.h>
#include <algorithm>

#include "ap_int.h"

extern int t_minimizer_sequence(
    const uint64_t* packed_sequence,
    const int n,
    const int s,
    const int w,
    ap_uint<64>* tab_hash
);

int main(int argc, char * argv[])
{
    printf("\n");
    printf("=============================================\n");
    printf(" Minimizer test   : %s\n", __FILE__);
    printf(" Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("=============================================\n\n");

    std::string i_file;
    std::string o_file = "result.txt";

    int s = 19;
    int w = 16;
    bool verbose = true;

 
    // -----------------------------
    // Parsing arguments
    // -----------------------------
    for (int i = 1; i < argc; i++) {
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
    // Lecture du fichier
    // -----------------------------
    std::string nucl;
    std::ifstream ii(i_file);

    if (ii.is_open()) {
        std::getline(ii, nucl);
        ii.close();
    }
    else {
        printf("(EE) Impossible to open input file (%s)\n", i_file.c_str());
        return -1;
    }

    const int n = nucl.size();

    if (n == 0) {
        printf("(EE) No data were loaded !\n");
        return -1;
    }

    const int n_smers     = n - s + 1;
    const int n_minzr_est = (2 * n_smers) / (w + 1);

    printf("------------ Parameters --------------------\n");
    printf(" Sequence length       : %d\n", n);
    printf(" s-mer size (s)        : %d\n", s);
    printf(" Window size (w)       : %d\n", w);
    printf(" Number of s-mers      : %d\n", n_smers);
    printf(" Estimated minimizers  : %d\n", n_minzr_est);
    printf("--------------------------------------------\n\n");

    // -------------------------
    // Allocation mémoire
    // -------------------------
    ap_uint<64>* tab_hash = new ap_uint<64>[n];

    // -------------------------
    // Appel accélérateur
    // -------------------------
    const int nHashs = t_minimizer_sequence(
        (const uint64_t*) nucl.c_str(),
        n,
        s,
        w,
        tab_hash
    );

    // -------------------------
    // Affichage des résultats
    // -------------------------
    int nPrint = 5;

    printf("\n✅ REAL number of minimizers found  : %d\n", nHashs);

    printf(">>> 5 FIRST HASH VALUES\n");
    for (int i = 0; i < std::min(nHashs, nPrint); ++i) {
        printf("  [%4d]  0x%016llX\n", i, tab_hash[i].to_uint64());
    }

    printf("\n>>> 5 LAST HASH VALUES\n");
    int start = (nHashs > nPrint) ? (nHashs - nPrint) : 0;
    for (int i = start; i < nHashs; ++i) {
        printf("  [%4d]  0x%016llX\n", i, tab_hash[i].to_uint64());
    }

    printf("\n");

    // -------------------------
    // Sauvegarde dans fichier
    // -------------------------
    FILE* ff = fopen(o_file.c_str(), "w");
    if (!ff) {
        printf("(EE) Impossible to open output file (%s)\n", o_file.c_str());
    }
    else {
        for (int i = 0; i < nHashs; ++i) {
            fprintf(ff, "[%4d] 0x%016llX\n", i, tab_hash[i].to_uint64());
        }
        fclose(ff);
        printf(" Results written in : %s\n", o_file.c_str());
    }

    delete[] tab_hash;

    printf("\n✅ DONE SUCCESSFULLY\n");

    return 0;
}
