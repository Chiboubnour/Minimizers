#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
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

    std::string i_file;
    std::string o_file = "result.txt";

    int s = 19;
    int w = 16;

    bool verbose = true;

    // ----------------------------------
    // Parse arguments
    // ----------------------------------
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

    if (i_file.empty())
    {
        printf("(EE) No input file provided (-i)\n");
        return EXIT_FAILURE;
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
        return EXIT_FAILURE;
    }

    const int n = nucl.size();
    if (n == 0)
    {
        printf("(EE) No data were loaded!\n");
        return EXIT_FAILURE;
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
    // 
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

    uint64_t total = nMinizrs.to_uint64();

    printf("\n✅ REAL number of minimizers found: %" PRIu64 "\n\n", total);

 
    if (verbose && total > 0)
    {
        const int show = 5;
        int first_count = (total < show) ? total : show;
        int last_count  = (total > show) ? show : 0;

        printf("-------- First %d minimizers --------\n", first_count);
        for (int i = 0; i < first_count; i++)
        {
            printf("s-mer [%6d] : hash 0x%016" PRIx64 "\n",
                   i,
                   tab_hash[i].to_uint64());
        }

        if (last_count > 0)
        {
            printf("-------- Last %d minimizers --------\n", last_count);
            for (uint64_t i = total - last_count; i < total; i++)
            {
                printf("s-mer [%6" PRIu64 "] : hash 0x%016" PRIx64 "\n",
                       i,
                       tab_hash[i].to_uint64());
            }
        }
    }

  
    FILE* ff = fopen(o_file.c_str(), "w");

    if (ff != nullptr)
    {
        fprintf(ff, "Number of minimizers found: %" PRIu64 "\n\n", total);

        for (uint64_t i = 0; i < total; ++i)
        {
            fprintf(ff,
                    "s-mer [%6" PRIu64 "] : hash 0x%016" PRIx64 "\n",
                    i,
                    tab_hash[i].to_uint64());
        }
        fclose(ff);

        printf("\n✅ Results successfully saved into: %s\n", o_file.c_str());
    }
    else
    {
        printf("(EE) Cannot open output file %s\n", o_file.c_str());
    }

    delete[] tab_hash;

    printf("\n✅ DONE\n");
    return 0;
}
