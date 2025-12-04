#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <inttypes.h>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cstdlib>

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
    int nTests = 1;
    bool verbose = true;

   
    for (int i = 1; i < argc; i++)
    {
        const std::string argvi = argv[i];

        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i")
        {
            i_file = argv[i + 1];
            i++;
        }
        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o")
        {
            o_file = argv[i + 1];
            i++;
        }
        else if (argvi == "--window" || argvi == "-w")
        {
            w = std::atoi(argv[i + 1]);
            i++;
        }
        else if (argvi == "--smer" || argvi == "-s")
        {
            s = std::atoi(argv[i + 1]);
            i++;
        }
        else if (argvi == "--t" || argvi == "-t")
        {
            nTests = std::atoi(argv[i + 1]);
            i++;
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

    if (i_file.empty())
    {
        printf("(EE) No input file provided (-i <file>)\n");
        exit(EXIT_FAILURE);
    }

    // -------------------------
    // Read input sequence
    // -------------------------
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
    const int n_smers     = n - s + 1;
    const int n_minzr_est = 2 * n_smers / (w + 1);

    printf("\n======== PARAMETERS ========\n");
    printf("Sequence size         : %d\n", n);
    printf("Number of s=%d-mers   : %d\n", s, n_smers);
    printf("Window size           : %d\n", w);
    printf("Estimated minimizers  : %d\n", n_minzr_est);
    printf("============================\n");


    ap_uint<64>* tab_hash = new ap_uint<64>[n];

    std::vector<double> durations_us;
    int nHashs = 0;

    for (int i = 0; i < nTests; ++i)
    {
        auto start = std::chrono::high_resolution_clock::now();

        nHashs = t_minimizer_sequence(
            (const uint64_t*) nucl.c_str(), 
            s,
            w,
            tab_hash
        );

        auto end = std::chrono::high_resolution_clock::now();

        double time_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        durations_us.push_back(time_us);

        if (verbose && i == 0)
        {
            printf("\n--- 5 first minimizers ---\n");
            for (int j = 0; j < nHashs && j < 5; ++j)
            {
                printf("s-mer [%2d] : hash 0x%016" PRIx64 "\n",
                       j, tab_hash[j].to_uint64());
            }

            printf("\n--- 5 last minimizers ---\n");
            for (int j = (nHashs > 5 ? nHashs - 5 : 0); j < nHashs; ++j)
            {
                printf("s-mer [%2d] : hash 0x%016" PRIx64 "\n",
                       j, tab_hash[j].to_uint64());
            }
        }

    }

    // -------------------------
    // Stats
    // -------------------------
    std::sort(durations_us.begin(), durations_us.end());

    double time_min    = durations_us.front() / 1e6;
    double time_max    = durations_us.back()  / 1e6;
    double time_median = durations_us[durations_us.size() / 2] / 1e6;
    double time_avg    = std::accumulate(durations_us.begin(),
                                         durations_us.end(), 0.0)
                          / (durations_us.size() * 1e6);

    double hashes_per_sec = nHashs / time_avg;
    double dHash           = hashes_per_sec / 1e6;

    const int mHash  = nHashs / (1024.0 * 1024.0);
    const int Mbytes = (nHashs * sizeof(uint64_t)) / (1024.0 * 1024.0);

    printf("\n✅ REAL number of minimizers found: %d\n", nHashs);

    std::cout << "\n#(II) Final results\n";
    std::cout << "#(II) - #of hash       : " << mHash   << " Mhash\n";
    std::cout << "#(II) - #of bytes      : " << Mbytes << " Mbytes\n";
    std::cout << "#(II) - Avg time       : " << time_avg    << " s\n";
    std::cout << "#(II) - Min time       : " << time_min    << " s\n";
    std::cout << "#(II) - Max time       : " << time_max    << " s\n";
    std::cout << "#(II) - Median time    : " << time_median << " s\n";
    std::cout << "#(II) - Throughput      : " << dHash << " M hash/s\n";

    printf("\n%4d %4d  %1.6f %1.6f %1.6f %1.6f %7.1f\n",
           Mbytes, mHash,
           time_avg, time_min,
           time_max, time_median,
           dHash);

    // -------------------------
    // Write output file
    // -------------------------
    FILE* ff = fopen(o_file.c_str(), "w");
    if (ff)
    {
        fprintf(ff, "5 first minimizers:\n");
        for (int i = 0; i < nHashs && i < 5; ++i)
        {
            fprintf(ff, "s-mer [%2d] : hash 0x%016" PRIx64 "\n",
                    i, tab_hash[i].to_uint64());
        }

        fprintf(ff, "\n5 last minimizers:\n");
        for (int i = (nHashs > 5 ? nHashs - 5 : 0); i < nHashs; ++i)
        {
            fprintf(ff, "s-mer [%2d] : hash 0x%016" PRIx64 "\n",
                    i, tab_hash[i].to_uint64());
        }

        fclose(ff);
        printf("\n✅ Results saved to file: %s\n", o_file.c_str());
    }
    else
    {
        printf("(EE) Impossible to open output file (%s)\n", o_file.c_str());
    }

    delete[] tab_hash;
    return EXIT_SUCCESS;
}
