#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring> 
#include <fstream> 
#include <inttypes.h>

#include "ap_int.h"

// Déclaration de la fonction HLS
extern int t_thread_minimizer(
    const uint64_t* packed_sequence,
    const int s,
    const int w,
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
    std::string o_file = "result_thr_v1.txt";

    int s = 19;
    int w = 16;
    bool verbose = true;

    // Affichage de la commande
    printf("# Run command:\n# ");
    for (int i = 0; i < argc; i++) printf("%s ", argv[i]);
    printf("\n");

    // Parsing des arguments
    for (int i = 1; i < argc; i++)
    {
        std::string argvi = argv[i];
        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i") {
            i_file = argv[++i];
        } else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o") {
            o_file = argv[++i];
        } else if (argvi == "--window" || argvi == "-w") {
            w = std::atoi(argv[++i]);
        } else if (argvi == "--smer" || argvi == "-s") {
            s = std::atoi(argv[++i]);
        } else if (argvi == "--no-verbose") {
            verbose = false;
        } else if (argvi == "-verbose") {
            verbose = true;
        } else {
            printf("(EE) Unknown argument: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    // Lecture de la séquence
    std::string nucl;
    std::ifstream ii(i_file);
    if(ii.is_open()){
        std::getline(ii, nucl);
        ii.close();
    } else {
        printf("(EE) Impossible to open input file (%s)\n", i_file.c_str());
        exit(EXIT_FAILURE);
    }

    const int n = nucl.size();
    if(n == 0){
        printf("(EE) No data were loaded !\n");
        exit(EXIT_FAILURE);
    }

    const int n_smers = n - s + 1;
    const int n_minzr = 2 * n_smers / (w + 1);

    printf("Taille de la sequence     : %d\n", n);
    printf("Nombre de s=%d-mers       : %d\n", s, n_smers);
    printf("Taille de la window.      : %d\n", w);
    printf("Nombre de (%d) minimizers : %d\n", s, n_minzr);

    // -------------------------
    // Pack ASCII sequence en uint64_t
    // -------------------------
    int n64 = (n + 7) / 8;
    std::vector<uint64_t> packed_seq(n64, 0);
    for(int i = 0; i < n; i++) {
        int word_idx = i / 8;
        int byte_idx = i % 8;
        packed_seq[word_idx] |= ((uint64_t)nucl[i] & 0xFF) << (8 * byte_idx);
    }

    // Allocation mémoire pour les hash
    ap_uint<64>* tab_hash = new ap_uint<64>[n]; // safe upper bound

    // Variable pour stocker le nombre de minimizers
    ap_uint<64> nMin = 0;

    // -------------------------
    // Appel de la fonction threadée
    // -------------------------
    const int nHashs = t_thread_minimizer(
        packed_seq.data(),
        s,
        w,
        tab_hash,
        &nMin
    );

    // Affichage
    if(verbose){
        for (int i = 0; i < nMin; ++i){
            printf("s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
        }
    }

    // Écriture dans le fichier
    FILE* ff = fopen(o_file.c_str(), "w");
    if(ff){
        for (int i = 0; i < nMin; ++i){
            fprintf(ff, "s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
        }
        fclose(ff);
    } else {
        printf("(EE) Impossible to open output file (%s)\n", o_file.c_str());
    }

    delete[] tab_hash;
    return 0;
}
