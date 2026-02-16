#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <ap_int.h>
#include "functions.hpp"


extern "C" {
void minimizer_v3(
    const uint64_t* seq,
    int s,
    int w,
    int n,
    uint64_t* tab_hash,
    uint64_t* nMinizrs
);
}

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

    // ---------- Parse arguments ----------
    for (int i = 1; i < argc; i++)
    {
        const std::string argvi = argv[i];

        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i")
            i_file = argv[++i];
        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o")
            o_file = argv[++i];
        else if (argvi == "--window" || argvi == "-w")
            w = std::atoi(argv[++i]);
        else if (argvi == "--smer" || argvi == "-s")
            s = std::atoi(argv[++i]);
        else if (argvi == "--no-verbose")
            verbose = false;
        else if (argvi == "--verbose")
            verbose = true;
        else {
            printf("(EE) Unknown argument: %s\n", argv[i]);
            exit(EXIT_FAILURE);
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
        nucl.clear();
        while (std::getline(ii, line)) {
            if (line.empty() || line[0]=='>') continue; // ignore FASTA headers
            nucl += line;
        }
        ii.close();
    } else {
        printf("(EE) Impossible to open input file (%s)\n", i_file.c_str());
        return EXIT_FAILURE;
    }

    const int n_seq = nucl.size();
    if (n_seq == 0) {
        printf("(EE) No data were loaded!\n");
        return EXIT_FAILURE;
    }

    const int n_smers = n_seq - s + 1;
    const int n_minzr_est = 2 * n_smers / (w + 1);

    printf("\n======== PARAMETERS ========\n");
    printf("Sequence size         : %d\n", n_seq);
    printf("Number of s=%d-mers   : %d\n", s, n_smers);
    printf("Window size           : %d\n", w);
    printf("Estimated minimizers  : %d\n", n_minzr_est);
    printf("============================\n");

    // ---------- Pack sequence into uint64_t words (8 bases per word) ----------
int n_words = (n_seq + 7) / 8;
uint64_t* seq64 = new uint64_t[n_words];
for (int i = 0; i < n_words; i++) {
    uint64_t word = 0;
    for (int j = 0; j < 8; j++) {
        int idx = i*8 + j;
        if (idx < n_seq) {
            uint8_t nuc2 = /*nucl_encode_v1(*/nucl[idx]/*)*/; // encode 2 bits
            word |= ((uint64_t)nuc2) << (8*j);
        }
    }
    seq64[i] = word;
}



    int tab_size = (n_smers + 7) / 8; // 8 s-mers per 512-bit word
    ap_uint<512>* tab_hash = new ap_uint<512>[tab_size];

    // ---------- Call HLS minimizer ----------
    uint64_t nMinizrs = 0;
    minimizer_v3(
        seq64,
        s,
        w,
        n_seq,
        reinterpret_cast<uint64_t*>(tab_hash),
        &nMinizrs
    );

    printf("\n✅ REAL number of minimizers found: %" PRIu64 "\n\n", nMinizrs);

    // ---------- Print all minimizers ----------
   if (verbose && nMinizrs > 0) {
    printf("-------- All minimizers (%" PRIu64 ") --------\n", nMinizrs);

    for (uint64_t i = 0; i < nMinizrs; i++) {
        uint64_t val = (uint64_t)tab_hash[i/8].range((i%8 + 1)*64 - 1, (i%8)*64);
        printf("Minimizer [%6" PRIu64 "] : hash 0x%016" PRIx64 "\n", i, val);
    }
    
   }

    // ---------- Save results to file ----------
    FILE* ff = fopen(o_file.c_str(), "w");
    if (ff != nullptr) {
        fprintf(ff, "Number of minimizers found: %" PRIu64 "\n\n", nMinizrs);
        for (uint64_t i = 0; i < nMinizrs; ++i) {
            uint64_t val = (uint64_t)tab_hash[i/8].range((i%8+1)*64-1, (i%8)*64);
            fprintf(ff, "Minimizer [%6" PRIu64 "] : hash 0x%016" PRIx64 "\n", i, val);
        }
        fclose(ff);
        printf("\n✅ Results successfully saved into: %s\n", o_file.c_str());
    } else {
        printf("(EE) Cannot open output file %s\n", o_file.c_str());
    }

    delete[] seq64;
    delete[] tab_hash;

    printf("\n✅ DONE\n");
    return 0;
}
