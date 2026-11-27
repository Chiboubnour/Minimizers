#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <cstdint>
#include <ap_int.h>

// ============================================================
// Fixed parameters (must match thr_minimizer_v3.cpp)
// ============================================================
constexpr int SMER   = 21;
constexpr int WINDOW = 11;

// ============================================================
// HLS kernel declaration
// Must match the implementation exactly
// ============================================================
void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs,
    ap_uint<64> n_bases
);

// ============================================================
// MAIN
// ============================================================
int main(int argc, char* argv[])
{
    printf("=====================================\n");
    printf(" Minimizer test : %s\n", __FILE__);
    printf(" Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("=====================================\n");

    std::string i_file;
    std::string o_file = "result.txt";

    bool verbose = true;

    // ---------- Parse arguments ----------
    for (int i = 1; i < argc; i++)
    {
        const std::string argvi = argv[i];

        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i") {
            if (++i >= argc) {
                printf("(EE) Missing input filename\n");
                return EXIT_FAILURE;
            }
            i_file = argv[i];
        }
        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o") {
            if (++i >= argc) {
                printf("(EE) Missing output filename\n");
                return EXIT_FAILURE;
            }
            o_file = argv[i];
        }
        else if (argvi == "--no-verbose") {
            verbose = false;
        }
        else if (argvi == "--verbose") {
            verbose = true;
        }
        else {
            printf("(EE) Unknown argument: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
    }

    if (i_file.empty()) {
        printf("(EE) No input file provided (-i)\n");
        return EXIT_FAILURE;
    }

    // ---------- Read input sequence ----------
    std::ifstream ii(i_file);
    std::string nucl;

    if (ii.is_open()) {
        std::string line;

        while (std::getline(ii, line)) {
            if (line.empty() || line[0] == '>')
                continue;

            nucl += line;
        }

        ii.close();
    }
    else {
        printf("(EE) Impossible to open input file (%s)\n",
               i_file.c_str());
        return EXIT_FAILURE;
    }

    const int n_seq = nucl.size();

    if (n_seq < SMER) {
        printf("(EE) Sequence too short!\n");
        return EXIT_FAILURE;
    }

    // ---------- Parameters ----------
    const int s = SMER;
    const int w = WINDOW;

    const int n_smers = n_seq - s + 1;
    const int n_minzr_est = 2 * n_smers / (w + 1);

    printf("\n======== PARAMETERS ========\n");
    printf("Sequence size         : %d\n", n_seq);
    printf("SMER                  : %d\n", s);
    printf("WINDOW                : %d\n", w);
    printf("Number of s-mers      : %d\n", n_smers);
    printf("Estimated minimizers  : %d\n", n_minzr_est);
    printf("============================\n");

    // ---------- Pack sequence into ap_uint<64> ----------
    // 8 bases per 64-bit word

    const int n_words = (n_seq + 7) / 8;

    ap_uint<64>* seq64 = new ap_uint<64>[n_words];

    for (int i = 0; i < n_words; i++) {

        ap_uint<64> word = 0;

        for (int j = 0; j < 8; j++) {

            int idx = i * 8 + j;

            if (idx < n_seq) {

                uint8_t nuc2 = static_cast<uint8_t>(nucl[idx]);

                word.range(8*j + 7, 8*j) = nuc2;
            }
        }

        seq64[i] = word;
    }

    // ---------- Allocate output ----------
    // One minimizer per ap_uint<64> element

    ap_uint<64>* tab_hash = new ap_uint<64>[n_smers];

    ap_uint<64> nMinizrs = 0;

    // ---------- Call HLS minimizer ----------
    printf("\nStarting minimizer kernel...\n");

    minimizer(
        seq64,
        tab_hash,
        &nMinizrs,
        ap_uint<64>(n_seq)
    );

    const uint64_t n_min = nMinizrs.to_uint64();

    printf("\nREAL number of minimizers found: %" PRIu64 "\n\n",
           n_min);

    // ---------- Print minimizers ----------
    if (verbose && n_min > 0) {

        printf("-------- All minimizers (%" PRIu64 ") --------\n",
               n_min);

        for (uint64_t i = 0; i < n_min; i++) {

            uint64_t val = tab_hash[i].to_uint64();

            printf("Minimizer [%6" PRIu64 "] : hash 0x%016" PRIx64 "\n",
                   i, val);
        }
    }

    // ---------- Save results ----------
    FILE* ff = fopen(o_file.c_str(), "w");

    if (ff != nullptr) {

        fprintf(ff,
                "Number of minimizers found: %" PRIu64 "\n\n",
                n_min);

        for (uint64_t i = 0; i < n_min; i++) {

            uint64_t val = tab_hash[i].to_uint64();

            fprintf(ff,
                    "Minimizer [%6" PRIu64 "] : hash 0x%016" PRIx64 "\n",
                    i, val);
        }

        fclose(ff);

        printf("\nResults successfully saved into: %s\n",
               o_file.c_str());
    }
    else {
        printf("(EE) Cannot open output file %s\n",
               o_file.c_str());
    }

    // ---------- Free memory ----------
    delete[] seq64;
    delete[] tab_hash;

    printf("\nDONE\n");

    return 0;
}
