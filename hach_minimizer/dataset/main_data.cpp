#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <cstdint>
#include <chrono>

#include <filesystem>
namespace fs = std::filesystem;


bool generate_minimizers_file(const fs::path& filepath, size_t size_mb) {
    const size_t bytes = size_mb * 1024 * 1024;
    const size_t count = bytes / sizeof(uint64_t);
    
    try {
        if (!filepath.parent_path().empty()) {
            fs::create_directories(filepath.parent_path());
        }
    } catch (const std::exception& e) {
        std::cerr << "Erreur lors de la création du répertoire: " << e.what() << std::endl;
        return false;
    }
    
    std::vector<uint64_t> buffer(count);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;

    std::cout << "Génération de " << count << " entiers ("
              << (bytes / (1024.0 * 1024.0)) << " MB) dans " << filepath << "..." << std::endl;

    // Remplir le buffer avec des nombres aléatoires
    #pragma omp parallel for
    for (size_t i = 0; i < count; ++i) {
        buffer[i] = dis(gen);
    }

    // Écrire les données dans le fichier
    std::ofstream outfile(filepath, std::ios::binary);
    if (!outfile) {
        std::cerr << "Erreur: Impossible d'écrire dans le fichier " << filepath.string() << std::endl;
        return false;
    }

    outfile.write(reinterpret_cast<const char*>(buffer.data()), bytes);
    if (!outfile) {
        std::cerr << "Erreur lors de l'écriture dans le fichier " << filepath.string() << std::endl;
        return false;
    }

    std::cout << "Fichier " << filepath.filename() << " généré avec succès (" 
              << (bytes / (1024.0 * 1024.0)) << " MB)" << std::endl;
    return true;
}

int main() {
    std::cout << "Génération de fichiers de minimizers aléatoires..." << std::endl;
    
    // Chemin de sortie pour les fichiers générés
    fs::path output_dir = "generated";
    
    // Créer le répertoire de sortie s'il n'existe pas
    try {
        fs::create_directories(output_dir);
    } catch (const std::exception& e) {
        std::cerr << "Erreur lors de la création du répertoire de sortie: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    bool success = true;
    success &= generate_minimizers_file(output_dir / "minimizers_32MB.bin",  32);
    success &= generate_minimizers_file(output_dir / "minimizers_64MB.bin",  64);
    success &= generate_minimizers_file(output_dir / "minimizers_128MB.bin", 128);
    success &= generate_minimizers_file(output_dir / "minimizers_256MB.bin", 256);
    
    if (success) {
        try {
            std::cout << "Génération terminée avec succès dans le répertoire: " 
                      << fs::absolute(output_dir).string() << std::endl;
        } catch (const std::exception& e) {
            std::cout << "Génération terminée avec succès dans le répertoire: " 
                      << output_dir.string() << std::endl;
            std::cerr << "Note: Impossible d'obtenir le chemin absolu: " << e.what() << std::endl;
        }
        return EXIT_SUCCESS;
    } else {
        std::cerr << "Des erreurs sont survenues lors de la génération des fichiers." << std::endl;
        return EXIT_FAILURE;
    }
}