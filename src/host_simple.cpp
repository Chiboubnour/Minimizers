#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cstring>

#include "xrt/xrt_bo.h"
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"

// Configuration
#define DATA_WIDTH 512
#define HASH_WIDTH 64
#define UNITS_PER_WORD (DATA_WIDTH / HASH_WIDTH)  // 8 unités de 64 bits

// Génération de données aléatoires
void generate_random_minimizers(std::vector<ap_uint<DATA_WIDTH>>& data, size_t num_words) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);
    
    data.resize(num_words);
    
    for (size_t i = 0; i < num_words; i++) {
        ap_uint<DATA_WIDTH> word = 0;
        
        // Génération de 8 minimizers de 64 bits chacun
        for (int j = 0; j < UNITS_PER_WORD; j++) {
            uint64_t random_value = dis(gen);
            word.range((j+1)*HASH_WIDTH-1, j*HASH_WIDTH) = random_value;
        }
        
        data[i] = word;
    }
}

// Affichage des résultats (premiers éléments)
void print_results(const std::vector<ap_uint<DATA_WIDTH>>& input, 
                   const std::vector<ap_uint<DATA_WIDTH>>& output, 
                   int num_to_show = 2) {
    
    std::cout << "\n=== Résultats (premiers " << num_to_show << " mots) ===\n";
    
    for (int i = 0; i < num_to_show && i < input.size(); i++) {
        std::cout << "\nMot " << i << ":\n";
        
        for (int j = 0; j < UNITS_PER_WORD; j++) {
            ap_uint<HASH_WIDTH> input_val = input[i].range((j+1)*HASH_WIDTH-1, j*HASH_WIDTH);
            ap_uint<HASH_WIDTH> output_val = output[i].range((j+1)*HASH_WIDTH-1, j*HASH_WIDTH);
            
            std::cout << "  Minimizer[" << j << "]: 0x" << std::hex << input_val.to_uint64() 
                      << " -> Hash: 0x" << output_val.to_uint64() << std::dec << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <xclbin_file> <device_id> <data_size_mb>\n";
        std::cout << "  data_size_mb: Taille des données en MB (ex: 512 pour 512MB)\n";
        return EXIT_FAILURE;
    }

    std::string binaryFile = argv[1];
    int device_index = std::stoi(argv[2]);
    int data_size_mb = std::stoi(argv[3]);
    
    // Calcul du nombre de mots de 512 bits
    size_t data_size_bytes = data_size_mb * 1024 * 1024;
    size_t num_words = data_size_bytes / (DATA_WIDTH / 8);  // 512 bits = 64 bytes
    
    std::cout << "=== Configuration ===\n";
    std::cout << "Taille des données: " << data_size_mb << " MB\n";
    std::cout << "Nombre de mots (512 bits): " << num_words << "\n";
    std::cout << "Minimizers par mot: " << UNITS_PER_WORD << "\n";
    std::cout << "Total minimizers: " << num_words * UNITS_PER_WORD << "\n\n";

    // Initialisation du device
    auto device = xrt::device(device_index);
    auto uuid = device.load_xclbin(binaryFile);
    auto krnl = xrt::kernel(device, uuid, "krnl_hash_simple");

    // Génération des données aléatoires
    std::cout << "Génération des données aléatoires...\n";
    auto start_gen = std::chrono::high_resolution_clock::now();
    
    std::vector<ap_uint<DATA_WIDTH>> input_data;
    generate_random_minimizers(input_data, num_words);
    
    auto end_gen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> gen_time = end_gen - start_gen;
    std::cout << "Génération terminée en " << gen_time.count() << " secondes\n\n";

    // Allocation des buffers HBM
    size_t buffer_size = num_words * sizeof(ap_uint<DATA_WIDTH>);
    
    auto bo_input = xrt::bo(device, buffer_size, XCL_BO_FLAGS_NONE, krnl.group_id(0));   // HBM[0]
    auto bo_output = xrt::bo(device, buffer_size, XCL_BO_FLAGS_NONE, krnl.group_id(1));  // HBM[1]

    auto input_map = bo_input.map<ap_uint<DATA_WIDTH>*>();
    auto output_map = bo_output.map<ap_uint<DATA_WIDTH>*>();

    // Copie des données vers HBM[0]
    std::cout << "Copie des données vers HBM[0]...\n";
    std::memcpy(input_map, input_data.data(), buffer_size);
    bo_input.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    // Exécution du kernel
    std::cout << "Exécution du kernel...\n";
    auto start_kernel = std::chrono::high_resolution_clock::now();
    
    auto run = krnl(bo_input, bo_output, static_cast<ap_uint<64>>(num_words));
    run.wait();
    
    auto end_kernel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> kernel_time = end_kernel - start_kernel;

    // Récupération des résultats depuis HBM[1]
    bo_output.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    std::vector<ap_uint<DATA_WIDTH>> output_data(num_words);
    std::memcpy(output_data.data(), output_map, buffer_size);

    // Affichage des résultats
    print_results(input_data, output_data, 3);
    
    // Statistiques de performance
    double total_hashes = num_words * UNITS_PER_WORD;
    double hashes_per_second = total_hashes / kernel_time.count();
    double throughput_mb_s = data_size_mb / kernel_time.count();
    
    std::cout << "\n=== Performance ===\n";
    std::cout << "Temps kernel: " << kernel_time.count() << " secondes\n";
    std::cout << "Hash calculés: " << total_hashes << "\n";
    std::cout << "Débit: " << hashes_per_second / 1e6 << " M hash/s\n";
    std::cout << "Throughput: " << throughput_mb_s << " MB/s\n";
    
    std::cout << "\nTraitement terminé avec succès!\n";
    return 0;
}
