#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <cstdint>
#include <string>
#include <cstdlib>
#include <numeric>   
#include <cstring>   
#include <algorithm> 
#include <getopt.h>
#include <filesystem>
namespace fs = std::filesystem;

// ta fonction HLS
void minimizer(
    const uint64_t* input_minimizers,
    uint64_t* output_hashes,
    unsigned int data_size_words
);

void read_uint64_t_file(
    const std::string& filename,
    std::vector<uint64_t>& i_buffer,
    const int nElements )
{
    fs::path filepath(filename);
    if (!fs::exists(filepath)) {
        std::cerr << "(EE) Unable to open file (" << filename << ")" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file) {
        std::cerr << "(EE) Unable to open file (" << filepath << ")" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    if (size % sizeof(uint64_t) != 0) {
        std::cerr << "(EE) The file size is not %8 (uint64_t)" << std::endl;
        exit(EXIT_FAILURE);
    }

    int nCases = size / sizeof(uint64_t);
    if (nElements != -1) nCases = std::min(nCases, nElements);

    i_buffer.resize(nCases);
    if (!file.read(reinterpret_cast<char*>(i_buffer.data()), nCases * sizeof(uint64_t))) {
        std::cerr << "(EE) An error happens during the data reading" << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char* argv[])
{
    int verbose_flag = 0;
    int help_flag = 0;
    int nElements = -1;
    int nTests = 64;
    std::string i_file = "", o_file = "";

    static struct option long_options[] = {
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"nhash", required_argument, 0, 'n'},
        {"ntest", required_argument, 0, 't'},
        {0,0,0,0}
    };

    int c, option_index = 0;
    while (true) {
        c = getopt_long(argc, argv, "i:o:n:t:vh", long_options, &option_index);
        if (c == -1) break;

        switch(c) {
            case 'i': i_file = optarg; break;
            case 'o': o_file = optarg; break;
            case 'n': nElements = 1024 * 1024 * std::atoi(optarg); break;
            case 't': nTests = std::atoi(optarg); break;
            case 'v': verbose_flag = true; break;
            case 'h': help_flag = true; break;
            default: abort();
        }
    }

    if (help_flag) {
        std::cout << "Usage: " << argv[0] << " -i <minimizer file> -o <output file>" << std::endl;
        return EXIT_FAILURE;
    }

    if (i_file.empty()) {
        std::cerr << "(EE) No input file provided..." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (o_file.empty()) {
        std::cerr << "(EE) No output file provided..." << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((nElements != -1) && (nElements < 1024*1024)) {
        std::cerr << "(EE) The number of hash to process is too small..." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (nTests < 1) {
        std::cerr << "(EE) The number of tests must be >= 1..." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (verbose_flag) std::cout << "#(II) Loading files from file (" << i_file << ")" << std::endl;

    std::vector<uint64_t> minimizers;
    read_uint64_t_file(i_file, minimizers, nElements);
    std::vector<uint64_t> hashes(minimizers.size());

    if (verbose_flag) std::cout << "#(II) Starting test procedure (" << nTests << " iter.)" << std::endl;

    std::vector<double> durations_us;
    durations_us.reserve(nTests);

    for (int i=0; i<nTests; ++i) {
        auto start = std::chrono::high_resolution_clock::now();

        // exécution de la fonction HLS
        minimizer(minimizers.data(), hashes.data(), minimizers.size());

        auto end = std::chrono::high_resolution_clock::now();
        durations_us.push_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());

        // affichage des premiers hashes pour debug
        if ((i==0) && verbose_flag) {
            std::cout << "#(II) Premiers hash générés :" << std::endl;
            for (size_t j=0; j<std::min<size_t>(5, hashes.size()); ++j)
                printf("#(II) - hash[%zu] = %16.16llX\n", j, hashes[j]);
        }

        // fake update pour éviter l'optimisation
        if (i<nTests-1)
            minimizers[rand() % minimizers.size()] = hashes[rand() % hashes.size()];
    }

    std::sort(durations_us.begin(), durations_us.end());
    double time_min    = durations_us.front() / 1e6;
    double time_max    = durations_us.back() / 1e6;
    double time_median = durations_us[durations_us.size()/2] / 1e6;
    double time_avg    = std::accumulate(durations_us.begin(), durations_us.end(), 0.0) / (durations_us.size() * 1e6);

    double total_hashes   = static_cast<double>(minimizers.size());
    double hashes_per_sec = total_hashes / time_avg;

    const int    mHash = minimizers.size() / 1e6;
    const double dHash = hashes_per_sec / 1e6;
    std::cout << "#(II) Final results" << std::endl;
    std::cout << "#(II) - #of hash      " <<  minimizers.size() / (1024.0*1024.0) << " Mhash" << std::endl;
    std::cout << "#(II) - #of bytes     " << (minimizers.size() * sizeof(uint64_t)) / (1024.0*1024.0) << " Mbytes" << std::endl;
    std::cout << "#(II) - Temps moyen:  " << time_avg << " s" << std::endl;
    std::cout << "#(II) - Temps min:    " << time_min << " s" << std::endl;
    std::cout << "#(II) - Temps max:    " << time_max << " s" << std::endl;
    std::cout << "#(II) - Temps médian: " << time_median << " s" << std::endl;
    std::cout << "#(II) - Débit:        " << dHash << " M hash/s" << std::endl;

    printf("%4d %1.6f %1.6f %1.6f %1.6f %7.1f\n", mHash, time_avg, time_min, time_max, time_median, dHash);

    return EXIT_SUCCESS;
}
