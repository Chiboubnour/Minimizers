#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <cstdint>
#include <string>
#include <algorithm> // pour std::sort
#include <numeric>   // pour std::accumulate


 void minimizer(
    const uint64_t* input_minimizers,
    uint64_t* output_hashes,
    unsigned int data_size_words
);

std::vector<uint64_t> read_minimizers_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    if (size % sizeof(uint64_t) != 0) {
        std::cerr << "Erreur: Taille de fichier invalide" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::vector<uint64_t> buffer(size / sizeof(uint64_t));
    if (!file.read(reinterpret_cast<char*>(buffer.data()), size)) {
        std::cerr << "Erreur lors de la lecture du fichier" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return buffer;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fichier_minimizers> [num_tests=10]" << std::endl;
        return EXIT_FAILURE;
    }

    std::string input_file = argv[1];
    int num_tests = (argc > 2) ? std::stoi(argv[2]) : 10;

    std::cout << "Chargement des minimizers depuis " << input_file << "..." << std::endl;
    auto minimizers = read_minimizers_file(input_file);
    std::vector<uint64_t> hashes(minimizers.size());

    std::cout << "Démarrage des tests avec " << num_tests << " itérations..." << std::endl;
    
    std::vector<double> durations_us;
    durations_us.reserve(num_tests);

    for (int i = 0; i < num_tests; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        
        minimizer(minimizers.data(), hashes.data(), minimizers.size());

        // Affichage au premier test uniquement
        if (i == 0) {
            std::cout << "\nPremiers minimizers en entrée :" << std::endl;
            for (size_t j = 0; j < std::min<size_t>(5, minimizers.size()); ++j) {
                std::cout << "minimizer[" << j << "] = " << minimizers[j] << std::endl;
            }

            std::cout << "\nPremiers hash générés :" << std::endl;
            for (size_t j = 0; j < std::min<size_t>(5, hashes.size()); ++j) {
                std::cout << "hash[" << j << "] = " << hashes[j] << std::endl;
            }
            std::cout << std::endl;
        }
        
        if (i < num_tests - 1) {
            minimizers[rand() % minimizers.size()] = hashes[rand() % hashes.size()];
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double time_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        durations_us.push_back(time_us);
        
    }

    std::sort(durations_us.begin(), durations_us.end());
    double time_min = durations_us.front() / 1e6;
    double time_max = durations_us.back() / 1e6;
    double time_median = durations_us[durations_us.size()/2] / 1e6;
    double time_avg = std::accumulate(durations_us.begin(), durations_us.end(), 0.0) / (durations_us.size() * 1e6);
    
    double total_hashes = minimizers.size();
    double hashes_per_sec = total_hashes / time_median;
    
    std::cout << "\nRésultats finaux:" << std::endl;
    std::cout << "  Taille des données: " << (minimizers.size() * sizeof(uint64_t)) / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  Temps min: " << time_min << " s" << std::endl;
    std::cout << "  Temps max: " << time_max << " s" << std::endl;
    std::cout << "  Temps médian: " << time_median << " s" << std::endl;
    std::cout << "  Temps moyen: " << time_avg << " s" << std::endl;
    std::cout << "  Débit: " << (hashes_per_sec / 1e6) << " M hash/s" << std::endl;
    
    return EXIT_SUCCESS;
}
