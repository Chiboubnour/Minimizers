#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

#define DATA_WIDTH 512
#define HASH_WIDTH 64
#define PARALLEL_UNITS 8
#define UNITS_PER_WORD (DATA_WIDTH / HASH_WIDTH)  // 8 unités de 64 bits dans 512 bits

inline ap_uint<64> bfc_hash_64(ap_uint<64> key, ap_uint<64> mask) {
    #pragma HLS INLINE
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

void parallel_hash_calc(
    ap_uint<DATA_WIDTH> input_data,
    ap_uint<DATA_WIDTH>& output_data
) {
    #pragma HLS INLINE
    
    ap_uint<HASH_WIDTH> hash_results[UNITS_PER_WORD];
    
    PARALLEL_HASH_LOOP: for (int i = 0; i < UNITS_PER_WORD; i++) {
        #pragma HLS UNROLL
        
        // Extraction de chaque unité de 64 bits
        ap_uint<HASH_WIDTH> input_unit = input_data.range(
            (i+1)*HASH_WIDTH-1, i*HASH_WIDTH
        );
        
        ap_uint<64> mask = ~0ULL;  // Masque complet pour 64 bits
        hash_results[i] = bfc_hash_64(input_unit, mask);
    }
    
    // Reconstruction du mot de 512 bits
    for (int i = 0; i < UNITS_PER_WORD; i++) {
        #pragma HLS UNROLL
        output_data.range((i+1)*HASH_WIDTH-1, i*HASH_WIDTH) = hash_results[i];
    }
}

extern "C" {
void krnl_hash_simple(
    const ap_uint<DATA_WIDTH>* input_minimizers,  
    ap_uint<DATA_WIDTH>* output_hashes,           
    ap_uint<64> data_size_words                   
) {
    #pragma HLS INTERFACE m_axi port=input_minimizers offset=slave bundle=gmem_in max_read_burst_length=16
    #pragma HLS INTERFACE m_axi port=output_hashes offset=slave bundle=gmem_out max_write_burst_length=16
    #pragma HLS INTERFACE s_axilite port=data_size_words
    #pragma HLS INTERFACE s_axilite port=return
    
    #pragma HLS DATAFLOW
    
    // Pipeline pour traiter les données par chunks
    const int BURST_SIZE = 16;  // Traitement par bursts de 16 mots
    
    MAIN_LOOP: for (ap_uint<64> burst_idx = 0; burst_idx < data_size_words; burst_idx += BURST_SIZE) {
        #pragma HLS PIPELINE II=1
        
        // Calcul de la taille du burst actuel
        ap_uint<64> remaining = data_size_words - burst_idx;
        ap_uint<64> current_burst = (remaining < BURST_SIZE) ? remaining : BURST_SIZE;
        
        // Traitement du burst
        BURST_LOOP: for (ap_uint<64> word_idx = 0; word_idx < current_burst; word_idx++) {
            #pragma HLS PIPELINE II=1
            #pragma HLS UNROLL factor=4
            
            ap_uint<64> addr = burst_idx + word_idx;
            
            // Lecture des données d'entrée
            ap_uint<DATA_WIDTH> input_data = input_minimizers[addr];
            
            // Calcul parallèle des hash
            ap_uint<DATA_WIDTH> output_data;
            parallel_hash_calc(input_data, output_data);
            
            // Écriture des résultats
            output_hashes[addr] = output_data;
        }
    }
}
}
