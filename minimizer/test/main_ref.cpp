#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring> 
#include <fstream> 
#include <inttypes.h>

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

    // Affichage commande
    printf("# Run command:\n# ");
    for (int i = 0; i < argc; i++) printf("%s ", argv[i]);
    printf("\n");

    // Arguments
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

    // Lecture du fichier
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

    ap_uint<64>* tab_hash = new ap_uint<64>[n]; // safe upper bound

    // Appel fonction avec 5 paramètres (conforme au .cpp)
    const int nHashs = t_minimizer_sequence(
        (const uint64_t*) nucl.c_str(),
        n, 
        s,
        w,
        tab_hash
    );

    if(verbose){
        for (int i = 0; i < nHashs; ++i){
            printf("s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
        }
    }

    FILE* ff = fopen(o_file.c_str(), "w");
    for (int i = 0; i < nHashs; ++i){
        fprintf(ff, "s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
    }
    fclose(ff);

    delete[] tab_hash;

    return 0;
}
