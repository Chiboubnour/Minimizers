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
    const int s,
    const int w,
    ap_uint<64>* tab_hash
);

int main(int argc, char * argv[])
{
    // -------------------------
    // Header
    // -------------------------
    printf("# \n");
    printf("# Minimizer test   : %s\n", __FILE__);
    printf("# Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("# \n");

    
    std::string i_file;
    std::string o_file = "result.txt";

    int s = 19; 
    int w = 16;  
    bool verbose = true;


    for (int i = 1; i < argc; i++)
    {
        const std::string argvi = argv[i];

        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i")
        {
            i_file = argv[i + 1];
            i += 1;
        }
        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o")
        {
            o_file = argv[i + 1];
            i += 1;
        }
        else if (argvi == "--window" || argvi == "-w")
        {
            w = std::atoi(argv[i + 1]);
            i += 1;
        }
        else if (argvi == "--smer" || argvi == "-s")
        {
            s = std::atoi(argv[i + 1]);
            i += 1;
        }
        else if (argvi == "--no-verbose")
        {
            verbose = false;
        }
        else if (argvi == "-verbose")
        {
            verbose = true;
        }
        else
        {
            printf("(EE) Unknown argument: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    std::string nucl;
    std::ifstream ii(i_file);
    if (ii.is_open())
    {
        std::getline(ii, nucl);
        ii.close();
    }
    else
    {
        printf("(EE) Impossible to open input file (%s)\n", i_file.c_str());
        exit(EXIT_FAILURE);
    }

    const int n = nucl.size();
    if (n == 0)
    {
        printf("(EE) No data were loaded!\n");
        exit(EXIT_FAILURE);
    }

    const int n_smers = n - s + 1;
    const int n_minzr_est = 2 * n_smers / (w + 1);

    printf("\n======== PARAMETERS ========\n");
    printf("Sequence size         : %d\n", n);
    printf("Number of s=%d-mers   : %d\n", s, n_smers);
    printf("Window size           : %d\n", w);
    printf("Estimated minimizers  : %d\n", n_minzr_est);
    printf("============================\n");

    ap_uint<64>* tab_hash = new ap_uint<64>[n];

    const int nHashs = t_minimizer_sequence(
        (const uint64_t*) nucl.c_str(),
        s, w,
        tab_hash
    );

    printf("\n✅ REAL number of minimizers found: %d\n", nHashs);

  
    printf("5 premiers hashes :\n");
    for (int i = 0; i < nHashs && i < 5; ++i)
    {
        printf("s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
    }

    printf("5 derniers hashes :\n");
    for (int i = (nHashs > 5 ? nHashs - 5 : 0); i < nHashs; ++i)
    {
        printf("s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
    }


    FILE* ff = fopen(o_file.c_str(), "w");
    if (ff)
    {
        fprintf(ff, "Number of minimizers found: %d\n\n", nHashs);

        fprintf(ff, "5 first hashes:\n");
        for (int i = 0; i < nHashs && i < 5; ++i)
        {
            fprintf(ff, "s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
        }

        fprintf(ff, "\n5 last hashes:\n");
        for (int i = (nHashs > 5 ? nHashs - 5 : 0); i < nHashs; ++i)
        {
            fprintf(ff, "s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
        }

        fclose(ff);
    }
    else
    {
        printf("(EE) Impossible to open output file (%s)\n", o_file.c_str());
    }

    delete[] tab_hash;

    return 0;
}
