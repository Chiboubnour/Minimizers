#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstring> 
#include <fstream> 
#include <inttypes.h>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "gz/reader/stream_gz_reader.hpp"
#include "ap_int.h"

extern int t_thread_minimizer(
    const uint64_t* packed_sequence,
    const int s,
    const int w,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
);

int main(int argc, char * argv[])
{
    printf("\n=============================================\n");
    printf(" Minimizer test   : %s\n", __FILE__);
    printf(" Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("=============================================\n\n");

    std::string i_file;                 
    std::string o_file = "result_thr_v1.txt";

    int s = 19;
    int w = 16;
    int ntests = 1; 
    bool verbose = true;
    uint64_t max_bases             = 0; 
    uint64_t total_bases_extracted = 0;

    // -----------------------------
    // Parsing des arguments
    // -----------------------------
    for (int i = 1; i < argc; i++)
    {
        std::string argvi = argv[i];

        if (argvi == "--input" || argvi == "--ifile" || argvi == "-i")
                i_file = argv[++i];

        else if (argvi == "--output" || argvi == "--ofile" || argvi == "-o") 
                o_file = argv[++i];

        else if (argvi == "--window" || argvi == "-w") 
                w = std::atoi(argv[++i]);

        else if (argvi == "--smer" || argvi == "-s") 
                s = std::atoi(argv[++i]);

        else if (argvi == "--ntests" || argvi == "-t") 
            ntests = std::atoi(argv[++i]);

        else if ((argvi == "--data") || (argvi == "-d")) {
            int mb = std::atoi(argv[++i]);
            max_bases = static_cast<uint64_t>(mb) * 1024ULL * 1024ULL;
        }

        else if (argvi == "--no-verbose") 
            verbose = false;

        else if (argvi == "-verbose") 
             verbose = true;

        else { printf("(EE) Unknown argument: %s\n", argv[i]); return -1; }
    }

    if (i_file.empty()) {
        printf("(EE) No input file provided. Use -i <file>\n");
        return -1;
    }

    constexpr int BLOCK_SIZE = 1024*1024; 
    std::vector<char> block_data;
    block_data.reserve(BLOCK_SIZE);

    std::vector<double> test_durations;
    int nHashs_total = 0;
    std::vector<ap_uint<64>> first_hashes, last_hashes;

    auto process_block = [&](std::vector<char> &block){
        if(block.empty()) return;
        int n = block.size();
        int n_smers = n - s + 1;
        if(n_smers <= 0) return;

        int n64 = (n + 7)/8;
        std::vector<uint64_t> packed(n64, 0ULL);
        for(int i=0;i<n;i++)
        {
            int word = i/8;
            int byte = i%8;
            packed[word] |= ((uint64_t)(uint8_t)block[i] << (8*byte));
        }

        ap_uint<64>* tab_hash = new ap_uint<64>[ n_smers ];
        ap_uint<64> nMin_hw = 0;

        auto start = std::chrono::high_resolution_clock::now();
        int nHashs = t_thread_minimizer(packed.data(), s, w, tab_hash, &nMin_hw);
        auto end = std::chrono::high_resolution_clock::now();

        double us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        test_durations.push_back(us / 1e6); // convert to seconds
        nHashs_total += nHashs;
        total_bases_extracted += n;

        for(int j=0; j<nHashs && first_hashes.size()<5; ++j)
            first_hashes.push_back(tab_hash[j]);
        for(int j=(nHashs>5 ? nHashs-5 :0); j<nHashs; ++j)
            last_hashes.push_back(tab_hash[j]);

        delete[] tab_hash;
        block.clear();
    };

    auto process_line = [&](std::string &line){
        if(line.empty() || line[0]=='>') return;
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        for(char c: line) block_data.push_back(c);
        if(block_data.size() >= BLOCK_SIZE) process_block(block_data);
    };

    auto read_file = [&](const std::string &file){
        block_data.clear();
        if(file.size()>=3 && file.substr(file.size()-3)==".gz"){
            stream_gz_reader gz(file);
            if(!gz.is_open()){ printf("(EE) Cannot open gz %s\n", file.c_str()); return; }
            std::vector<char> buf(64*1024);
            std::string partial;
            partial.reserve(1024);
            while(!gz.is_eof())
            {
                int r = gz.read(buf.data(),1,(int)buf.size());
                if(r<=0) break;
                partial.append(buf.data(), r);
                size_t pos = 0;
                while(true)
                {
                    size_t nl = partial.find('\n', pos);
                    if(nl==std::string::npos){
                        partial = partial.substr(pos);
                        break;
                    }
                    std::string line = partial.substr(pos, nl-pos);
                    if(!line.empty() && line.back()=='\r') line.pop_back();
                    process_line(line);
                    pos = nl+1;
                }
            }
            if(!partial.empty()) process_line(partial);
            gz.close();
        } else {
            std::ifstream ifs(file);
            if(!ifs.is_open()){ printf("(EE) Cannot open %s\n", file.c_str()); return; }
            std::string line;
            while(std::getline(ifs,line))
            {
                if(!line.empty() && line.back()=='\r') line.pop_back();
                process_line(line);
            }
            if(!block_data.empty()) process_block(block_data);
            ifs.close();
        }
    };

    for(int t=0; t<ntests; ++t){
        read_file(i_file);
    }

    // ==== Affichage des statistiques ====
    if(!test_durations.empty())
    {
        std::sort(test_durations.begin(), test_durations.end());
        double time_min    = test_durations.front();
        double time_max    = test_durations.back();
        double time_median = test_durations[test_durations.size()/2];
        double time_avg    = std::accumulate(test_durations.begin(), test_durations.end(),0.0)/test_durations.size();
        double hashes_per_sec = nHashs_total / time_avg;
        double dHash           = hashes_per_sec / 1e6;

        std::cout << "\n# Total bases read   : " << total_bases_extracted << "\n";
        std::cout << "# Total minimizers   : " << nHashs_total << "\n";
        std::cout << "# Min time           : " << time_min << " s\n";
        std::cout << "# Max time           : " << time_max << " s\n";
        std::cout << "# Median time        : " << time_median << " s\n";
        std::cout << "# Average time       : " << time_avg << " s\n";
        std::cout << "# Throughput         : " << dHash << " M hash/s\n";
    }

    std::cout << "\n# First 5 minimizers:\n";
    for (size_t i=0;i<first_hashes.size() && i<5; ++i)
        printf("  [%zu] 0x%016" PRIx64 "\n", i, first_hashes[i].to_uint64());

    std::cout << "\n# Last 5 minimizers:\n";
    for (size_t i=(last_hashes.size()>5?last_hashes.size()-5:0); i<last_hashes.size(); ++i)
        printf("  [%zu] 0x%016" PRIx64 "\n", i, last_hashes[i].to_uint64());

    // Écriture dans fichier
    FILE* ff = fopen(o_file.c_str(),"w");
    if(ff){
        for(int i=0;i<nHashs_total;++i) fprintf(ff,"0x%016llX\n",first_hashes[i].to_uint64());
        fclose(ff);
        printf("\n✅ Results written to %s\n", o_file.c_str());
    }

    printf("\n✅ DONE SUCCESSFULLY\n");
    return 0;
}
