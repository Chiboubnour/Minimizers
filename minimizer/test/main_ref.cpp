#include <algorithm>
#include <chrono>
#include <cmath>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>


template <int s, int w>
int t_minimizer_sequence(
    const char* sequence,
    int         n_bases,
    uint64_t*   tab_hash
);


int run_minimizer(
    const char* sequence,
    int         n_bases,
    int         s,
    int         w,
    uint64_t*   tab_hash
)
{
    // Add cases here if you want to benchmark other (s,w)
    if (s == 19 && w == 16)
        return t_minimizer_sequence<19,16>(sequence, n_bases, tab_hash);

    std::cerr << "(EE) Unsupported parameters: s=" << s
              << " w=" << w << '\n';
    std::exit(EXIT_FAILURE);
}

// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[])
{
    std::printf("# \n");
    std::printf("# Minimizer benchmark : %s\n", __FILE__);
    std::printf("# Compilation date    : %s %s\n", __DATE__, __TIME__);
    std::printf("# \n");

    // -------------------------
    // Parameters (defaults)
    // -------------------------
    std::string i_file;
    std::string o_file = "result.txt";

    int  s        = 19;
    int  w        = 16;
    int  ntests   = 1;
    bool verbose  = true;

    uint64_t max_bases             = 0; // 0 = no limit
    uint64_t total_bases_extracted = 0;

    // -------------------------
    // Argument parsing
    // -------------------------
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--input") || (arg == "--ifile") || (arg == "-i")) {
            i_file = argv[++i];
        }
        else if ((arg == "--output") || (arg == "--ofile") || (arg == "-o")) {
            o_file = argv[++i];
        }
        else if ((arg == "--window") || (arg == "-w")) {
            w = std::atoi(argv[++i]);
        }
        else if ((arg == "--smer") || (arg == "-s")) {
            s = std::atoi(argv[++i]);
        }
        else if ((arg == "--ntests") || (arg == "-t")) {
            ntests = std::atoi(argv[++i]);
        }
        else if (arg == "--no-verbose") {
            verbose = false;
        }
        else if (arg == "--verbose") {
            verbose = true;
        }
        else if ((arg == "--data") || (arg == "-d")) {
            int mb = std::atoi(argv[++i]);
            max_bases = static_cast<uint64_t>(mb) * 1024ULL * 1024ULL;
        }
        else {
            std::printf("(EE) Unknown argument: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
    }

    if (i_file.empty()) {
        std::printf("(EE) No input file provided (-i)\n");
        return EXIT_FAILURE;
    }

    if (ntests < 1) {
        std::printf("(EE) ntests must be >= 1\n");
        return EXIT_FAILURE;
    }

 
    std::vector<char>      sequence;
    std::vector<uint64_t>  first_hashes;
    std::vector<double>    test_durations;

    int nHashs_total = 0;

    // -------------------------
    // Read & prepare input
    // -------------------------
    auto read_and_prepare_data = [&]() -> bool {
        sequence.clear();
        total_bases_extracted = 0;

        auto process_line = [&](std::string& line) {
            if (line.empty() || line[0] == '>') return;

            line.erase(
                std::remove_if(line.begin(), line.end(), ::isspace),
                line.end()
            );

            for (char c : line) {
                if (max_bases > 0 && total_bases_extracted >= max_bases)
                    return;

                sequence.push_back(c);
                ++total_bases_extracted;
            }
        };

        std::ifstream ifs(i_file);
        if (!ifs.is_open()) {
            std::printf("(EE) Cannot open %s\n", i_file.c_str());
            return false;
        }

        std::string line;
        while (std::getline(ifs, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            process_line(line);
            if (max_bases > 0 && total_bases_extracted >= max_bases) break;
        }
        ifs.close();

        return true;
    };

    if (!read_and_prepare_data()) return EXIT_FAILURE;

    if (sequence.size() < static_cast<size_t>(s)) {
        std::cerr << "(EE) Sequence too small (input too small?). Exiting.\n";
        return EXIT_FAILURE;
    }

    // -------------------------
    // Benchmark loop
    // -------------------------
    test_durations.clear();
    test_durations.reserve(ntests);

    int reference_hash_count = -1;

    const int n_bases = static_cast<int>(sequence.size());

    for (int t = 0; t < ntests; ++t) {
        uint64_t* tab_hash = new uint64_t[n_bases];

        auto start = std::chrono::high_resolution_clock::now();
        int nHashs = run_minimizer(
            sequence.data(),
            n_bases,
            s,
            w,
            tab_hash
        );
        auto end = std::chrono::high_resolution_clock::now();

        double us = static_cast<double>(
            std::chrono::duration_cast<std::chrono::microseconds>(
                end - start
            ).count()
        );

        double time_s = us / 1e6;
        test_durations.push_back(time_s);

        if (t == 0) {
            first_hashes.assign(tab_hash, tab_hash + nHashs);

            FILE* ff = fopen(o_file.c_str(), "w");
            if (ff != nullptr) {
                for (int i = 0; i < nHashs; ++i) {
                    fprintf(ff, "Minimizer [%6d] : hash 0x%016" PRIx64 "\n",
                            i, first_hashes[i]);
                }
                fclose(ff);
                printf("\n✅ Results successfully saved into: %s\n", o_file.c_str());
            } else {
                printf("(EE) Cannot open output file %s\n", o_file.c_str());
            }

            reference_hash_count = nHashs;
            nHashs_total         = nHashs;
        }
        else if (nHashs != reference_hash_count) {
            std::cerr << "Warning: inconsistent number of hashes in test "
                      << t << ": expected " << reference_hash_count
                      << ", got " << nHashs << '\n';
        }

        delete[] tab_hash;
    }

    // -------------------------
    // Statistics
    // -------------------------
    if (!test_durations.empty()) {
        std::vector<double> stats;
        if (ntests > 1)
            stats.assign(test_durations.begin() + 1, test_durations.end());
        else
            stats = test_durations;
    
            if (!stats.empty()) {
                std::sort(stats.begin(), stats.end());
            
                double min_time = stats.front();
                double max_time = stats.back();
                double sum_time = std::accumulate(stats.begin(), stats.end(), 0.0);
                double avg_time = sum_time / static_cast<double>(stats.size());
            
                double var = 0.0;
                for (double v : stats) var += (v - avg_time) * (v - avg_time);
                var /= static_cast<double>(stats.size());
            
                double avg_hashes_per_sec = (avg_time > 0.0)
                    ? static_cast<double>(nHashs_total) / avg_time
                    : 0.0;
            
                double avg_bases_per_sec = (avg_time > 0.0)
                    ? static_cast<double>(total_bases_extracted) / avg_time
                    : 0.0;
            
                // ---- Detailed info ----
                printf("\n=== Performance Analysis ===\n\n");
                printf("SizeMB  nbases       t_moy       t_min        t_max       Mbases/s    NbMinz\n");
                printf("%6.1f  %9" PRIu64 "  %10.6f  %10.6f  %10.6f  %10.3f  %10" PRIu64 "\n",
                       static_cast<double>(total_bases_extracted)/(1024.0*1024.0),  // SizeMB
                       total_bases_extracted,                                        // nbases
                       avg_time,                                                      // t_moy
                       min_time,                                                      // t_min
                       max_time,                                                      // t_max
                       avg_bases_per_sec/1e6,                                         // Mbases/s
                       nHashs_total);                                                 // NbMinz
            
            
                    }
    }
    

    std::printf("\n Test completed successfully (output file: %s)\n",
                o_file.c_str());
    
    return EXIT_SUCCESS;
}
    


