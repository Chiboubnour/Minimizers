#if defined(__ARM_NEON__) || defined(__ARM_NEON)
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdint>
#include <arm_neon.h>
#include <arm_neon.h>
#include <stdint.h>
#include <stddef.h>
//
#include "../../tools/TimeMeasure.hpp"

// LUT (identique à avant)
static uint8_t shuffle_table[65536][16];
static uint8_t popcount_table[65536];

static void init_shuffle_table() {
    for (uint32_t m = 0; m < 65536; m++) {
        uint8_t c = 0;
        for (int i = 0; i < 16; i++) {
            if (m & (1 << i)) shuffle_table[m][c++] = i;
        }
        for (int i = c; i < 16; i++) shuffle_table[m][i] = 0;
        popcount_table[m] = c;
    }
}

static inline uint16_t movemask(uint8x16_t v) {
    uint8_t tmp[16];
    vst1q_u8(tmp, vshrq_n_u8(v, 7));
    uint16_t m = 0;
    for (int i = 0; i < 16; i++) m |= (tmp[i] & 1) << i;
    return m;
}

static inline uint8x16_t valid_mask(uint8x16_t v) {
    uint8x16_t A = vceqq_u8(v, vdupq_n_u8('A'));
    uint8x16_t C = vceqq_u8(v, vdupq_n_u8('C'));
    uint8x16_t G = vceqq_u8(v, vdupq_n_u8('G'));
    uint8x16_t T = vceqq_u8(v, vdupq_n_u8('T'));
    uint8x16_t U = vceqq_u8(v, vdupq_n_u8('U'));
    uint8x16_t N = vceqq_u8(v, vdupq_n_u8('N'));

    return vorrq_u8(vorrq_u8(vorrq_u8(A, C), vorrq_u8(G, T)), vorrq_u8(U, N));
}

// ===============================
// VERSION CORRIGÉE
// ===============================

void fasta_neon_ultimate(
    const char *in,
    size_t size,
    char *out,
    size_t *out_size
) {
    size_t i = 0, pos = 0;

    uint8_t in_header = 0;
    uint8_t writing = 0;

    for (; i + 15 < size; i += 16) {

        uint8x16_t v = vld1q_u8((const uint8_t*)(in + i));

        uint8x16_t is_gt = vceqq_u8(v, vdupq_n_u8('>'));
        uint8x16_t is_nl = vceqq_u8(v, vdupq_n_u8('\n'));

        // uppercase
        uint8x16_t is_lower = vandq_u8(
            vcgtq_u8(v, vdupq_n_u8('a' - 1)),
            vcltq_u8(v, vdupq_n_u8('z' + 1))
        );
        uint8x16_t upper = vsubq_u8(v, vandq_u8(is_lower, vdupq_n_u8(32)));

        uint8x16_t valid = valid_mask(upper);

        uint16_t gt_mask = movemask(is_gt);
        uint16_t nl_mask = movemask(is_nl);

        // ===== header mask =====
        uint16_t header_mask = 0;
        uint8_t state = in_header;

        for (int j = 0; j < 16; j++) {
            if (gt_mask & (1 << j)) state = 1;
            header_mask |= (state << j);
            if (nl_mask & (1 << j)) state = 0;
        }
        in_header = state;

        // ===== keep mask =====
        uint8x16_t keep_vec = vandq_u8(valid, vmvnq_u8(is_nl));
        uint16_t keep_mask = movemask(keep_vec);
        keep_mask &= ~header_mask;

        // ===== 🚨 NEWLINE FIX =====
        // on insère newline EXACTEMENT à la position du '>'
        // et AVANT les données suivantes

        for (int j = 0; j < 16; j++) {
            if (gt_mask & (1 << j)) {
                if (writing) {
                    out[pos++] = '\n';
                    writing = 0;
                }
            }

            if (keep_mask & (1 << j)) {
                // accès direct au byte SIMD
                out[pos++] = ((uint8_t*)&upper)[j];
                writing = 1;
            }
        }
    }

    // ===== reste scalaire =====
    for (; i < size; i++) {
        char c = in[i];

        if (c == '>') {
            if (writing) {
                out[pos++] = '\n';
                writing = 0;
            }
            in_header = 1;
        }
        else if (in_header) {
            if (c == '\n') in_header = 0;
        }
        else {
            if (c!=' ' && c!='\n' && c!='\t' && c!='\r') {
                char u = (c >= 'a' && c <= 'z') ? c - 32 : c;

                if (u=='A'||u=='C'||u=='G'||u=='T'||u=='U'||u=='N') {
                    out[pos++] = u;
                    writing = 1;
                }
            }
        }
    }

    if (writing) out[pos++] = '\n';

    *out_size = pos;
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
void fasta_neon_ultimate_parser(
    const std::string& file_i,
    const std::string& file_o
){
    init_shuffle_table();
    //
    //
    //
    FILE* f_in = fopen(file_i.c_str(), "rb");
    if( f_in == NULL ){
        printf("(EE) Error while opening input file [%s]\n", file_i.c_str());
        exit(EXIT_FAILURE);
    }
    fseek(f_in, 0, SEEK_END);
    size_t size = ftell(f_in);
    rewind(f_in);

    //
    //
    //
    FILE* f_ou = fopen(file_o.c_str(), "wb");
    if( f_ou == NULL ){
        printf("(EE) Error while opening output file [%s]\n", file_o.c_str());
        exit(EXIT_FAILURE);
    }

    //
    //
    //
    TimeMeasure c_time;

    c_time.start();
    char* buffer_i = (char*)malloc(size);
    int64_t rBytes = fread(buffer_i, 1, size, f_in);
    fclose(f_in);
    c_time.stop();

    std::cout << "Temps : " << c_time.ms() << " ms and debit : " << c_time.MBps(rBytes) << " Mbps" << std::endl;

    //
    // buffer sortie (taille max = entrée + nb séquences)
    //
    char* buffer_o = (char*)malloc(size + 1024);

    //
    //
    //
    c_time.start();
    size_t out_size;
    fasta_neon_ultimate(buffer_i, size, buffer_o, &out_size);
    c_time.stop();

    std::cout << "Temps : " << c_time.ms() << " ms and debit : " << c_time.MBps(rBytes) << " Mbps" << std::endl;
    
    //
    //
    //
    c_time.start();
    fwrite(buffer_o, 1, out_size, f_ou);
    fclose(f_ou);
    c_time.stop();

    std::cout << "Temps : " << c_time.ms() << " ms and debit : " << c_time.MBps(rBytes) << " Mbps" << std::endl;

    //
    //
    //
    free(buffer_i);
    free(buffer_o);
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
#endif