#if defined(__ARM_NEON__) || defined(__ARM_NEON)
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdint>
#include <stdint.h>
#include <stddef.h>
#include <arm_neon.h>
//
#include "../../tools/TimeMeasure.hpp"

// ==========================
// Validation nucléotides
// ==========================

static inline uint8x16_t valid_mask(const uint8x16_t v) {
    const uint8x16_t A = vceqq_u8(v, vdupq_n_u8('A'));
    const uint8x16_t C = vceqq_u8(v, vdupq_n_u8('C'));
    const uint8x16_t G = vceqq_u8(v, vdupq_n_u8('G'));
    const uint8x16_t T = vceqq_u8(v, vdupq_n_u8('T'));
    const uint8x16_t U = vceqq_u8(v, vdupq_n_u8('U'));
    const uint8x16_t N = vceqq_u8(v, vdupq_n_u8('N'));
    return vorrq_u8(vorrq_u8(vorrq_u8(A, C), vorrq_u8(G, T)), vorrq_u8(U, N));
}

// ==========================
// movemask NEON
// ==========================

inline uint16_t neon_movemask(const uint8x16_t input)
{
    const uint16x8_t high_bits = vreinterpretq_u16_u8 (vshrq_n_u8 (input, 7));
    const uint32x4_t paired16  = vreinterpretq_u32_u16(vsraq_n_u16(high_bits, high_bits, 7));
    const uint64x2_t paired32  = vreinterpretq_u64_u32(vsraq_n_u32(paired16, paired16, 14));
    const uint8x16_t paired64  = vreinterpretq_u8_u64 (vsraq_n_u64(paired32, paired32, 28));
    return vgetq_lane_u8(paired64, 0) | ((int)vgetq_lane_u8(paired64, 8) << 8);
}

// ==========================
// FONCTION PRINCIPALE
// ==========================

void fasta_neon(
    const char *in,
    size_t size,
    char *out,
    size_t *out_size
) {
    size_t i = 0, pos = 0;

    int in_header = 0;
    int writing = 0;

    for (; i + 15 < size; i += 16) {

        uint8x16_t v = vld1q_u8((const uint8_t*)(in + i));

        // Détection caractères
        uint8x16_t is_gt = vceqq_u8(v, vdupq_n_u8('>'));
        uint8x16_t is_nl = vceqq_u8(v, vdupq_n_u8('\n'));

        // Uppercase ASCII (branchless)
        uint8x16_t is_lower = vandq_u8(
            vcgtq_u8(v, vdupq_n_u8('a' - 1)),
            vcltq_u8(v, vdupq_n_u8('z' + 1))
        );
        uint8x16_t upper = vsubq_u8(v, vandq_u8(is_lower, vdupq_n_u8(32)));

        // Validité ADN/ARN
        uint8x16_t valid = valid_mask(upper);

        // Export vers tableau (pour logique header correcte)
        uint8_t buf[16];
        uint8_t m_valid[16];
        uint8_t m_gt[16];
        uint8_t m_nl[16];

        vst1q_u8(buf, upper);
        vst1q_u8(m_valid, valid);
        vst1q_u8(m_gt, is_gt);
        vst1q_u8(m_nl, is_nl);

        for (int j = 0; j < 16; j++) {

            if (m_gt[j]) {
                if (writing) {
                    out[pos++] = '\n';
                    writing = 0;
                }
                in_header = 1;
                continue;
            }

            if (in_header) {
                if (m_nl[j]) {
                    in_header = 0;
                }
                continue;
            }

            if (m_nl[j]) continue;

            if (m_valid[j]) {
                out[pos++] = buf[j];
                writing = 1;
            }
        }
    }

    // Reste scalaire
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

    if (writing) {
        out[pos++] = '\n';
    }

    *out_size = pos;
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
void fasta_neon_parser(
    const std::string& file_i,
    const std::string& file_o
){
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
    fasta_neon(buffer_i, size, buffer_o, &out_size);
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