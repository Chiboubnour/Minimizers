#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring> 
#include <fstream> 
#include <inttypes.h>

#include "ap_int.h"


 void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n_bases,
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

    std::string i_file ;
    std::string o_file = "result.txt";

    int s = 19;
    int w = 16;

    bool  verbose = true;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //
    printf("# Run command:\n# ");
    for (int i = 0; i < argc; i += 1)
    {
        printf("%s ", argv[i]);
    }
    printf("\n");
    //
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    for (int i = 1; i < argc; i++)
    {
        const std::string argvi = argv[i];
        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i" )
        {
            i_file = argv[i + 1];
            i += 1;
        }
        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o" )
        {
            o_file = argv[i + 1];
            i += 1;
        }
        else if (argvi == "--window" || argvi == "-w" )
        {
            w = std::atoi(argv[i + 1]);
            i += 1;
        }
        else if (argvi == "--smer" || argvi == "-s" )
        {
            s = std::atoi(argv[i + 1]);
            i += 1;
        }
        else if (std::string(argv[i]) == "--no-verbose")
        {
            verbose = false;
        }
        else if (std::string(argv[i]) == "-verbose")
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
    std::ifstream ii( i_file );
    if( ii.is_open() == true ){
        std::getline(ii, nucl);
        ii.close();
    }else{
        printf("(EE) Impossible to open input file (%s)\n", i_file.c_str() );
        exit(EXIT_FAILURE);
    }

    const int n = nucl.size();
    if( n == 0) {
        printf("(EE) No data were loaded !\n");
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

    //
    // Compression des nucleotides par lot de 64 bits
    //
    ap_uint<64>* packed_sequence = new ap_uint<64>[ (n+7) / 8 ];
    for(int i = 0; i < n; i += 8)
    {
        ap_uint<64> data = 0;
        for(int j = 0; j < 8; j += 1)
        {
            if( (i+j) < n )
                data = data | (((uint64_t)nucl[i+j]) << (8 *j));
        }
        packed_sequence[i/8] = data;
    }

    ap_uint<64>* tab_hash = new ap_uint<64>[ n ];
    ap_uint<64> nMinizrs;
    
    minimizer(
        packed_sequence,
        n,
        tab_hash,
        &nMinizrs
    );


        printf("\n✅ REAL number of minimizers found : %llu\n\n",
        (unsigned long long) nMinizrs.to_uint64());
        uint64_t total = nMinizrs.to_uint64();

        int first_count = (total < 5) ? total : 5;
        int last_count  = (total < 5) ? 0 : 5;

        printf("----- First %d minimizers -----\n", first_count);
        for (int i = 0; i < first_count; i++)
        {
            printf("s-mer [%6d] : hash 0x%016llX\n",
                   i,
                   (unsigned long long) tab_hash[i].to_uint64());
        }

        
        printf("----- Last %d minimizers -----\n", last_count);
        for (uint64_t i = total - 5; i < total; i++)
        {
                printf("s-mer [%6llu] : hash 0x%016llX\n",
                       (unsigned long long)i,
                       (unsigned long long)tab_hash[i].to_uint64());
        }
        
    

   
    FILE* ff = fopen(o_file.c_str(), "w");
    for (int i = 0; i < nMinizrs; ++i)
    {
        fprintf(ff, "s-mer [%2d] : hash 0x%16.16llX\n", i, tab_hash[i].to_uint64());
    }
    fclose( ff );

    delete[] tab_hash;

    return 0;
}

