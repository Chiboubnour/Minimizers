#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <ap_int.h>

void t_thread_minimizer_v2(
    const uint64_t* seq,
    int s,
    int w,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
);

int main(int argc, char * argv[])
{
    printf("=====================================\n");
    printf(" Minimizer test   : %s\n", __FILE__);
    printf(" Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("=====================================\n");

    std::string i_file ;
    std::string o_file = "result.txt";

    int s = 19;
    int w = 16;

    bool verbose = true;


    for (int i = 1; i < argc; i++)
    {
        const std::string argvi = argv[i];

        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i")
        {
            i_file = argv[++i];
        }
        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o")
        {
            o_file = argv[++i];
        }
        else if (argvi == "--window" || argvi == "-w")
        {
            w = std::atoi(argv[++i]);
        }
        else if (argvi == "--smer" || argvi == "-s")
        {
            s = std::atoi(argv[++i]);
        }
        else if (argvi == "--no-verbose")
        {
            verbose = false;
        }
        else if (argvi == "--verbose")
        {
            verbose = true;
        }
        else
        {
            printf("(EE) Unknown argument: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    // ----------------------------------
    // Read input sequence
    // ----------------------------------
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

    // ----------------------------------
    // Allocation large enough
    // ----------------------------------
    ap_uint<64>* tab_hash = new ap_uint<64>[n_smers];
    ap_uint<64> nMinizrs = 0;

    // ----------------------------------
    // HLS call
    // ----------------------------------
    t_thread_minimizer_v2(
        (const uint64_t*)nucl.c_str(),
        s,
        w,
        tab_hash,
        &nMinizrs
    );

    printf("\n✅ REAL number of minimizers found : %llu\n\n",
           (unsigned long long) nMinizrs.to_uint64());

  
    if (verbose && nMinizrs > 0)
    {
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

        if (last_count > 0)
        {
            printf("----- Last %d minimizers -----\n", last_count);
            for (uint64_t i = total - 5; i < total; i++)
            {
                printf("s-mer [%6llu] : hash 0x%016llX\n",
                       (unsigned long long)i,
                       (unsigned long long)tab_hash[i].to_uint64());
            }
        }
    }

  
    FILE* ff = fopen(o_file.c_str(), "w");

    if (ff != nullptr)
    {
        for (uint64_t i = 0; i < nMinizrs; ++i)
        {
            fprintf(ff, "s-mer [%6llu] : hash 0x%016llX\n",
                    (unsigned long long)i,
                    (unsigned long long)tab_hash[i].to_uint64());
        }
        fclose(ff);
    }
    else
    {
        printf("(EE) Cannot open output file %s\n", o_file.c_str());
    }

    delete[] tab_hash;

    printf("\n✅ DONE\n");
    return 0;
}
