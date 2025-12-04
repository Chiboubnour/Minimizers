// test/main_hls_v1.cpp
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <inttypes.h>

#include "ap_int.h"

// Le top HLS exporté par ton fichier HLS
extern "C" void minimizer_hls_top(
    const ap_uint<64>* packed_sequence,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
);

int main(int argc, char * argv[])
{
    // header
    printf("=============================================\n");
    printf("   HLS MINIMIZER TEST — %s\n", __FILE__);
    printf("   Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("=============================================\n\n");

    std::string i_file;
    std::string o_file = "result_hls.txt";

    int s = 19;
    int w = 16;
    bool verbose = true;

    // show run
    printf("# Run command:\n# ");
    for (int i = 0; i < argc; ++i) printf("%s ", argv[i]);
    printf("\n\n");

    // parse args (simple, safe)
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "--input" || a == "-i" || a == "--ifile") && i+1 < argc) { i_file = argv[++i]; }
        else if ((a == "--output" || a == "-o" || a == "--ofile") && i+1 < argc) { o_file = argv[++i]; }
        else if ((a == "--window" || a == "-w") && i+1 < argc) { w = std::atoi(argv[++i]); }
        else if ((a == "--smer" || a == "-s") && i+1 < argc) { s = std::atoi(argv[++i]); }
        else if (a == "--no-verbose") { verbose = false; }
        else if (a == "-verbose") { verbose = true; }
        else {
            // ignore unknown options but print
            printf("(WW) Unknown/ignored arg: %s\n", argv[i]);
        }
    }

    if (i_file.empty()) {
        printf("(EE) No input file provided. Use -i <file>\n");
        return -1;
    }

    // read first line
    std::string nucl;
    std::ifstream fi(i_file);
    if (!fi.is_open()) {
        printf("(EE) Cannot open input file: %s\n", i_file.c_str());
        return -1;
    }
    std::getline(fi, nucl);
    fi.close();

    const int n = static_cast<int>(nucl.size());
    if (n <= 0) {
        printf("(EE) Input sequence is empty.\n");
        return -1;
    }

    const int n_smers = n - s + 1;
    const int n_minzr_est = 2 * n_smers / (w + 1);

    printf("Sequence length : %d\n", n);
    printf("s-mer size (s)  : %d\n", s);
    printf("window (w)      : %d\n", w);
    printf("n s-mers        : %d\n", n_smers);
    printf("estimated minzs : %d\n", n_minzr_est);

    // ---------------------------------------------------------------------
    // HLS reader expects up to 64 words (64*8 = 512 chars). Build array of 64.
    // ---------------------------------------------------------------------
    const int MAX_WORDS = 64;
    const int n_words = (n + 7) / 8;
    ap_uint<64>* packed_sequence = new ap_uint<64>[MAX_WORDS];
    // init to 0 (important)
    for (int i = 0; i < MAX_WORDS; ++i) packed_sequence[i] = 0;

    // pack LSB-first (exactly as reference)
    for (int wi = 0; wi < n_words && wi < MAX_WORDS; ++wi) {
        ap_uint<64> w64 = 0;
        for (int j = 0; j < 8; ++j) {
            int idx = wi * 8 + j;
            if (idx < n) {
                uint64_t c = static_cast<uint8_t>(nucl[idx]);
                w64 |= ( (uint64_t)c << (8 * j) );
            }
        }
        packed_sequence[wi] = w64;
    }
    if (n_words > MAX_WORDS) {
        printf("WARNING: sequence truncated to %d chars (HLS reader limit = %d words)\n",
               MAX_WORDS*8, MAX_WORDS);
    }

    // allocate output buffer (safe upper bound = n)
    const int max_out = (n > 0) ? n : 1;
    ap_uint<64>* tab_hash = new ap_uint<64>[max_out];
    // zero-init to avoid reading uninitialized memory later
    for (int i = 0; i < max_out; ++i) tab_hash[i] = 0;

    ap_uint<64> nMinizrs = 0;

    // Safety check: ensure pointers non-null
    if (!packed_sequence || !tab_hash) {
        printf("(EE) Allocation failed\n");
        delete[] tab_hash;
        delete[] packed_sequence;
        return -1;
    }

    // Try/catch not useful for HLS C wrapper, but we print debug info before call.
    printf("Calling minimizer_hls_top(...)\n");
    minimizer_hls_top(packed_sequence, tab_hash, &nMinizrs);

    uint64_t found = nMinizrs.to_uint64();
    printf("HLS returned %llu minimizers\n", (unsigned long long)found);

    // Print a nice summary + first/last 5
    const int NPRINT = 5;
    printf("\n=== First %d minimizers ===\n", std::min<int>(NPRINT, (int)found));
    for (uint64_t i = 0; i < std::min<uint64_t>(NPRINT, found); ++i) {
        printf("s-mer [%2llu] : hash 0x%16.16llX\n",
               (unsigned long long)i,
               (unsigned long long)tab_hash[i].to_uint64());
    }

    printf("\n=== Last %d minimizers ===\n", std::min<int>(NPRINT, (int)found));
    uint64_t start = (found > NPRINT) ? (found - NPRINT) : 0;
    for (uint64_t i = start; i < found; ++i) {
        printf("s-mer [%2llu] : hash 0x%16.16llX\n",
               (unsigned long long)i,
               (unsigned long long)tab_hash[i].to_uint64());
    }

    // write file
    FILE* ff = fopen(o_file.c_str(), "w");
    if (ff) {
        for (uint64_t i = 0; i < found; ++i) {
            fprintf(ff, "s-mer [%2llu] : hash 0x%16.16llX\n",
                    (unsigned long long)i,
                    (unsigned long long)tab_hash[i].to_uint64());
        }
        fclose(ff);
        printf("\nResults saved to: %s\n", o_file.c_str());
    } else {
        printf("(EE) Cannot open output file %s\n", o_file.c_str());
    }

    delete[] tab_hash;
    delete[] packed_sequence;
    printf("\nDone.\n");
    return 0;
}
