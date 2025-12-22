// ============================================================
// Minimizer benchmark & validation tool (CPU reference)
// Compatible with t_minimizer_sequence<s,w>() provided earlier
// - Command-line interface
// - FASTA / FASTA.gz input
// - Statistics & performance metrics
// ============================================================

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

#include "gz/reader/stream_gz_reader.hpp"

// ------------------------------------------------------------
// Forward declaration of the reference minimizer (template)
// ------------------------------------------------------------

template <int s, int w>
int t_minimizer_sequence(
    const char* sequence,
    int         n_bases,
    uint64_t*   tab_hash
);

// ------------------------------------------------------------
// Runtime wrapper (dispatch on (s,w))
// ------------------------------------------------------------

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

    constexpr int BLOCK_SIZE = 1024 * 1024; // bases per block

    // -------------------------
    // Block structure
    // -------------------------
    struct BlockData {
        std::vector<char> sequence;
        int               n_bases = 0;
    };

    std::vector<BlockData> blocks;
    std::vector<uint64_t>  first_hashes;
    std::vector<uint64_t>  last_hashes;
    std::vector<double>    test_durations;

    int nHashs_total = 0;

    // -------------------------
    // Read & prepare input
    // -------------------------
    auto read_and_prepare_data = [&]() -> bool {
        std::vector<char> block_data;
        block_data.reserve(BLOCK_SIZE);

        total_bases_extracted = 0;
        blocks.clear();

        auto process_block = [&](std::vector<char>& block) {
            if (block.size() < static_cast<size_t>(s)) {
                block.clear();
                return;
            }
            BlockData bd;
            bd.sequence = std::move(block);
            bd.n_bases  = static_cast<int>(bd.sequence.size());
            blocks.push_back(std::move(bd));
            block.clear();
        };

        auto process_line = [&](std::string& line) {
            if (line.empty() || line[0] == '>') return;

            line.erase(
                std::remove_if(line.begin(), line.end(), ::isspace),
                line.end()
            );

            for (char c : line) {
                if (max_bases > 0 && total_bases_extracted >= max_bases) {
                    process_block(block_data);
                    return;
                }

                block_data.push_back(c);
                ++total_bases_extracted;

                if (static_cast<int>(block_data.size()) >= BLOCK_SIZE) {
                    process_block(block_data);
                }
            }
        };

        // --- Open input (supports .gz) ---
        if (i_file.size() >= 3 && i_file.substr(i_file.size() - 3) == ".gz") {
            stream_gz_reader gz(i_file);
            if (!gz.is_open()) {
                std::printf("(EE) Cannot open gz %s\n", i_file.c_str());
                return false;
            }

            std::vector<char> buf(64 * 1024);
            std::string       partial;
            partial.reserve(1024);

            while (!gz.is_eof()) {
                int r = gz.read(buf.data(), 1, static_cast<int>(buf.size()));
                if (r <= 0) break;

                partial.append(buf.data(), r);
                size_t pos = 0;

                while (true) {
                    size_t nl = partial.find('\n', pos);
                    if (nl == std::string::npos) {
                        partial = partial.substr(pos);
                        break;
                    }

                    std::string line = partial.substr(pos, nl - pos);
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    process_line(line);
                    pos = nl + 1;
                }
            }

            if (!partial.empty()) process_line(partial);
            gz.close();
        }
        else {
            std::ifstream ifs(i_file);
            if (!ifs.is_open()) {
                std::printf("(EE) Cannot open %s\n", i_file.c_str());
                return false;
            }

            std::string line;
            while (std::getline(ifs, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                process_line(line);
            }
            ifs.close();
        }

        if (!block_data.empty()) process_block(block_data);
        return true;
    };

    if (!read_and_prepare_data()) return EXIT_FAILURE;

    if (blocks.empty()) {
        std::cerr << "(EE) No blocks prepared (input too small?). Exiting.\n";
        return EXIT_FAILURE;
    }

    // -------------------------
    // Benchmark loop
    // -------------------------
    test_durations.clear();
    test_durations.reserve(ntests);

    int reference_hash_count = -1;

    for (int t = 0; t < ntests; ++t) {
        int current_test_hashes = 0;
        std::vector<double> durations_us;

        for (const auto& block : blocks) {
            uint64_t* tab_hash = new uint64_t[block.n_bases];

            auto start = std::chrono::high_resolution_clock::now();
            int nHashs = run_minimizer(
                block.sequence.data(),
                block.n_bases,
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

            durations_us.push_back(us);
            current_test_hashes += nHashs;

            if (t == 0) {
                for (int j = 0;
                     j < nHashs && static_cast<int>(first_hashes.size()) < 5;
                     ++j) {
                    first_hashes.push_back(tab_hash[j]);
                }
                for (int j = std::max(0, nHashs - 5); j < nHashs; ++j) {
                    last_hashes.push_back(tab_hash[j]);
                }
            }

            delete[] tab_hash;
        }

        double total_us = std::accumulate(
            durations_us.begin(), durations_us.end(), 0.0
        );
        double time_s = total_us / 1e6;
        test_durations.push_back(time_s);

        if (t == 0) {
            reference_hash_count = current_test_hashes;
            nHashs_total         = reference_hash_count;
        }
        else if (current_test_hashes != reference_hash_count) {
            std::cerr << "Warning: inconsistent number of hashes in test "
                      << t << ": expected " << reference_hash_count
                      << ", got " << current_test_hashes << '\n';
        }
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
                printf("Configuration:\n");
                printf("  Input file:          %s\n", i_file.c_str());
                printf("  Total bases:         %" PRIu64 "\n", total_bases_extracted);
                printf("  Window size (w):     %d\n", w);
                printf("  S-mer size (s):      %d\n", s);
                printf("  Number of tests:     %d\n\n", ntests);
            
                printf("Results (on %zu runs):\n", stats.size());
                printf("  Min time:            %.6f s\n", min_time);
                printf("  Max time:            %.6f s\n", max_time);
                printf("  Avg time:            %.6f s\n", avg_time);
            
                printf("Performance:\n");
                printf("  Number of minimizers: %" PRIu64 "\n", nHashs_total);
                printf("  Hash rate (avg):      %.3f Mhash/s\n", avg_hashes_per_sec / 1e6);
                printf("  Base rate (avg):      %.3f Mbases/s\n\n", avg_bases_per_sec / 1e6);
            
                // ---- One-line summary ----
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
    
    int first_count = std::min<int>(5, nHashs_total);
    int last_count  = std::min<int>(5, nHashs_total);
    uint64_t total  = nHashs_total;
    
    printf("\n----- First %zu minimizers -----\n", first_hashes.size());
    for (size_t i = 0; i < first_hashes.size(); ++i) {
        printf("s-mer [%6zu] : hash 0x%016llX\n",
               i,
               (unsigned long long)first_hashes[i]);
    }
    
    printf("\n----- Last %zu minimizers -----\n", last_hashes.size());
    for (size_t i = 0; i < last_hashes.size(); ++i) {
        printf("s-mer [%6zu] : hash 0x%016llX\n",
               i,
               (unsigned long long)last_hashes[i]);
    }
    
    std::printf("\n✅ Test completed successfully (output file: %s)\n",
                o_file.c_str());
    
    return EXIT_SUCCESS;
}
    