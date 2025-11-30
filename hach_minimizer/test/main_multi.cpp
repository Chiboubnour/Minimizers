#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <cstdint>
#include <string>
#include <cstdlib>
#include <omp.h>
#include <numeric> // Pour std::accumulate
#include <cstring> // Pour strerror
#include <numeric>   // pour std::accumulate
#include <algorithm> // pour std::sort


#include <filesystem>
namespace fs = std::filesystem;


// Déclaration de la fonction de hachage
extern "C" void minimizer_cpu_multi(const uint64_t* input_minimizers, uint64_t* output_hashes, unsigned int data_size_words);

// Fonction pour lire un fichier de minimizers
std::vector<uint64_t> read_minimizers_file(const std::string& filename) {
    fs::path filepath(filename);
    if (!fs::exists(filepath)) {
        std::cerr << "Erreur: Le fichier " << filename << " n'existe pas" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier " << filepath << " : " 
                  << strerror(errno) << std::endl;
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
        std::cerr << "Usage: " << argv[0] << " <fichier_minimizers> [num_tests=64] [num_threads=4]" << std::endl;
        return EXIT_FAILURE;
    }

    std::string input_file = argv[1];
    int num_tests = (argc > 2) ? std::stoi(argv[2]) : 64;   // 64 tests par défaut
    int num_threads = (argc > 3) ? std::stoi(argv[3]) : 4;  // 4 threads par défaut

    // Configurer le nombre de threads OpenMP
    omp_set_num_threads(num_threads);
    std::cout << "Utilisation de " << num_threads << " threads" << std::endl;

    // Résolution du chemin du fichier
    fs::path input_path(input_file);
    if (input_path.is_relative()) {
        if (fs::exists(input_path)) {
            input_path = fs::absolute(input_path);
        } else if (fs::exists("../dataset/" + input_file)) {
            input_path = fs::absolute(fs::path("../dataset") / input_file);
        } else if (fs::exists("../generated/" + input_file)) {
            input_path = fs::absolute(fs::path("../generated") / input_file);
        } else {
            std::cerr << "Erreur : le fichier " << input_file << " n'a pas été trouvé." << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "Chargement des minimizers depuis " << input_path << "..." << std::endl;
    auto minimizers = read_minimizers_file(input_path.string());
    std::vector<uint64_t> hashes(minimizers.size());

    std::cout << "Démarrage des tests avec " << num_tests << " itérations..." << std::endl;

    std::vector<double> durations_us;
    durations_us.reserve(num_tests);

    for (int i = 0; i < num_tests; ++i) {
        auto start = std::chrono::high_resolution_clock::now();

        // Appel de la fonction de hachage multi-thread
        minimizer_cpu_multi(minimizers.data(), hashes.data(), minimizers.size());

        if (i == 0) {
            std::cout << "\nPremiers hash générés :" << std::endl;
            for (size_t j = 0; j < std::min<size_t>(5, hashes.size()); ++j) {
                std::cout << "hash[" << j << "] = " << hashes[j] << std::endl;
            }
            std::cout << std::endl;
        }
    
        auto end = std::chrono::high_resolution_clock::now();
        double time_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        durations_us.push_back(time_us);


    }

    // Calcul des statistiques
    std::sort(durations_us.begin(), durations_us.end());
    double time_min = durations_us.front() / 1e6;
    double time_max = durations_us.back() / 1e6;
    double time_median = durations_us[durations_us.size() / 2] / 1e6;
    double time_avg = std::accumulate(durations_us.begin(), durations_us.end(), 0.0) / (durations_us.size() * 1e6);

    double total_hashes = static_cast<double>(minimizers.size());
    double hashes_per_sec = total_hashes / time_median;

    std::cout << "\n===== Résultats finaux =====" << std::endl;
    std::cout << "  Taille des données: " << (minimizers.size() * sizeof(uint64_t)) / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  Temps min: " << time_min << " s" << std::endl;
    std::cout << "  Temps max: " << time_max << " s" << std::endl;
    std::cout << "  Temps médian: " << time_median << " s" << std::endl;
    std::cout << "  Temps moyen: " << time_avg << " s" << std::endl;
    std::cout << "  Débit: " << (hashes_per_sec / 1e6) << " M hash/s" << std::endl;

    return EXIT_SUCCESS;
}
