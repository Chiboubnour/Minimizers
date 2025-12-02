#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <cstdint>
#include <filesystem>
namespace fs = std::filesystem;

bool generate_minimizers_file(const fs::path& filepath, size_t size_mb, size_t block_mb = 64) {
    const size_t total_bytes = size_mb * 1024 * 1024;
    const size_t block_bytes = block_mb * 1024 * 1024;
    const size_t uint64_per_block = block_bytes / sizeof(uint64_t);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;

    try {
        if (!filepath.parent_path().empty())
            fs::create_directories(filepath.parent_path());
    } catch (const std::exception& e) {
        std::cerr << "Erreur lors de la création du répertoire: " << e.what() << std::endl;
        return false;
    }

    std::ofstream outfile(filepath, std::ios::binary);
    if (!outfile) {
        std::cerr << "Erreur: impossible d'ouvrir " << filepath << " pour écriture" << std::endl;
        return false;
    }

    size_t bytes_written = 0;
    std::vector<uint64_t> buffer(uint64_per_block);

    while (bytes_written < total_bytes) {
        size_t bytes_to_write = std::min(block_bytes, total_bytes - bytes_written);

        for (size_t i = 0; i < bytes_to_write / sizeof(uint64_t); ++i) {
            buffer[i] = dis(gen);
        }

        outfile.write(reinterpret_cast<const char*>(buffer.data()), bytes_to_write);
        if (!outfile) {
            std::cerr << "Erreur lors de l'écriture dans " << filepath << std::endl;
            return false;
        }

        bytes_written += bytes_to_write;
        std::cout << "\rProgression: " << (bytes_written * 100 / total_bytes) << "%   " << std::flush;
    }

    std::cout << "\nFichier " << filepath.filename() << " généré avec succès (" 
              << size_mb << " MB)" << std::endl;
    return true;
}

int main() {
    fs::path output_file = "generated/minimizers_4GB.bin";

    if (generate_minimizers_file(output_file, 4096)) {
        std::cout << "Génération terminée avec succès dans : " << fs::absolute(output_file) << std::endl;
        return EXIT_SUCCESS;
    } else {
        std::cerr << "Erreur lors de la génération du fichier." << std::endl;
        return EXIT_FAILURE;
    }
}
