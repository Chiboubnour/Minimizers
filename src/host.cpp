#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>

#include "xrt/xrt_bo.h"
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"


std::vector<std::string> read_fasta_all(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier FASTA : " + filename);
    }

    std::vector<std::string> sequences;
    std::string line;
    std::string current_seq;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                sequences.push_back(std::move(current_seq));
                current_seq.clear();
            }
        } else {
            for (char ch : line) {
                char nuc = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                if (nuc == 'A' || nuc == 'C' || nuc == 'G' || nuc == 'T') {
                    current_seq.push_back(nuc);
                }
            }
        }
    }
    if (!current_seq.empty()) sequences.push_back(std::move(current_seq));
    return sequences;
}

inline uint64_t nucl_encode(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 0;
    }
}

// === Kernel runner ===
double run_krnl_reuse(xrt::device& device, xrt::kernel& krnl,
    xrt::bo& bo_seq, xrt::bo& bo_hash, xrt::bo& bo_nMinizrs,
    uint64_t* seq_map_ptr, uint64_t* hash_map_ptr, uint64_t* nmin_map_ptr,
    const std::vector<uint64_t>& packed_seq, size_t n)
{
    size_t input_size_bytes  = packed_seq.size() * sizeof(uint64_t);
    size_t n_smers           = (n >= 28) ? (n - 27) : 0; // S=28
    size_t output_size_bytes = n_smers * sizeof(uint64_t);

    // copy packed sequence into mapped pointer (only needed words)
    std::memcpy(seq_map_ptr, packed_seq.data(), input_size_bytes);

    // zero only the relevant portion of hash output
    std::memset(hash_map_ptr, 0, output_size_bytes);

    // clear count
    *nmin_map_ptr = 0;

    // sync to device (only input and control BO if needed)
    bo_seq.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_nMinizrs.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // launch kernel and measure kernel runtime only
    auto t0 = std::chrono::high_resolution_clock::now();
    auto run = krnl(bo_seq, static_cast<int>(n), bo_hash, bo_nMinizrs);
    run.wait();
    auto t1 = std::chrono::high_resolution_clock::now();

    bo_hash.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    bo_nMinizrs.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

    std::chrono::duration<double> kernel_time = t1 - t0;
    return kernel_time.count();
}

// === Main ===
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <xclbin_file> <device_id> <fasta_file>\n";
        return EXIT_FAILURE;
    }

    std::string binaryFile = argv[1];
    int device_index = std::stoi(argv[2]);
    std::string fastaFile = argv[3];

    std::cout << "Device: " << device_index << "\nLoading xclbin: " << binaryFile << std::endl;

    auto device = xrt::device(device_index);
    auto uuid = device.load_xclbin(binaryFile);
    auto krnl = xrt::kernel(device, uuid, "krnl_minimizer");

    std::vector<std::string> sequences;
    try {
        sequences = read_fasta_all(fastaFile);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    if (sequences.empty()) {
        std::cerr << "Aucune sequence trouvée dans le FASTA." << std::endl;
        return EXIT_FAILURE;
    }

    size_t max_n = 0;
    for (auto &s : sequences) if (s.size() > max_n) max_n = s.size();
    size_t max_n_words = (max_n + 7) / 8;
    size_t max_n_smers  = (max_n >= 28) ? (max_n - 27) : 0;
    size_t max_input_bytes  = max_n_words * sizeof(uint64_t);
    size_t max_output_bytes = max_n_smers * sizeof(uint64_t);

    std::cout << "Nombre de sequences lues : " << sequences.size() << "\n";
    std::cout << "Taille maximale sequence : " << max_n << " bases\n";

    
    int arg_index_seq = 0;
    int arg_index_hash = 2;
    int arg_index_nmin = 3;

    auto bo_seq = xrt::bo(device, max_input_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(arg_index_seq));
    auto bo_hash = xrt::bo(device, max_output_bytes, XCL_BO_FLAGS_NONE, krnl.group_id(arg_index_hash));
    auto bo_nMinizrs = xrt::bo(device, sizeof(uint64_t), XCL_BO_FLAGS_NONE, krnl.group_id(arg_index_nmin));

    auto seq_map_ptr = bo_seq.map<uint64_t*>();
    auto hash_map_ptr = bo_hash.map<uint64_t*>();
    auto nmin_map_ptr = bo_nMinizrs.map<uint64_t*>();

    size_t seq_idx = 0;
    for (const auto &sequence : sequences) {
        size_t n = sequence.size();
        size_t n_words = (n + 7) / 8;
        size_t n_smers = (n >= 28) ? (n - 27) : 0;
        size_t input_bytes = n_words * sizeof(uint64_t);
        size_t output_bytes = n_smers * sizeof(uint64_t);

        std::cout << "\n--- Sequence " << seq_idx << " : " << n << " bases ---\n";

        std::vector<uint64_t> packed_seq(n_words, 0);
        for (size_t i = 0; i < n; ++i) {
            size_t word_idx = i / 8;
            size_t shift = 8 * (i % 8);
            packed_seq[word_idx] |= (static_cast<uint64_t>(static_cast<uint8_t>(sequence[i])) << shift);
        }

        double ktime = run_krnl_reuse(device, krnl,
                                      bo_seq, bo_hash, bo_nMinizrs,
                                      seq_map_ptr, hash_map_ptr, nmin_map_ptr,
                                      packed_seq, n);

        uint64_t nmin_found = *nmin_map_ptr;

        std::cout << "Minimizers trouvés : " << nmin_found << "\n";
        std::cout << "Temps kernel (s) : " << ktime << "\n";

        ++seq_idx;
    }

    std::cout << "\nTraitement termine.\n";
    return 0;
}
