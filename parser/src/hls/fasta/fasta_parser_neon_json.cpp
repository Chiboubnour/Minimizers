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
#include "../../tools/stats.hpp"
//
//
//
//
//
// ================= LUT =================

static uint8_t shuffle_table[65536][16];
static uint8_t popcount_table[65536];

void init_shuffle_table() {
    for (uint32_t m = 0; m < 65536; m++) {
        uint8_t c = 0;
        for (int i = 0; i < 16; i++) {
            if (m & (1 << i)) {
                shuffle_table[m][c++] = i;
            }
        }
        for (int i = c; i < 16; i++) {
            shuffle_table[m][i] = 0;
        }
        popcount_table[m] = c;
    }
}

// ================= MOVEMASK =================

inline int movemask(const uint8x16_t input)
{
    uint16x8_t high_bits = vreinterpretq_u16_u8 (vshrq_n_u8 (input, 7));
    uint32x4_t paired16  = vreinterpretq_u32_u16(vsraq_n_u16(high_bits, high_bits, 7));
    uint64x2_t paired32  = vreinterpretq_u64_u32(vsraq_n_u32(paired16, paired16, 14));
    uint8x16_t paired64  = vreinterpretq_u8_u64 (vsraq_n_u64(paired32, paired32, 28));
    return vgetq_lane_u8(paired64, 0) | ((int)vgetq_lane_u8(paired64, 8) << 8);
}

// ================= VALID =================

static inline uint8x16_t valid_mask(uint8x16_t v) {
    uint8x16_t A = vceqq_u8(v, vdupq_n_u8('A'));
    uint8x16_t C = vceqq_u8(v, vdupq_n_u8('C'));
    uint8x16_t G = vceqq_u8(v, vdupq_n_u8('G'));
    uint8x16_t T = vceqq_u8(v, vdupq_n_u8('T'));
    uint8x16_t U = vceqq_u8(v, vdupq_n_u8('U'));
    uint8x16_t N = vceqq_u8(v, vdupq_n_u8('N'));

    return vorrq_u8(
        vorrq_u8(vorrq_u8(A, C), vorrq_u8(G, T)),
        vorrq_u8(U, N)
    );
}

// ================= MAIN =================

void fasta_neon_json(
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

        if (keep_mask == 0 && gt_mask == 0) continue;

        // ===============================
        // segmentation intra-bloc
        // ===============================

        int start = 0;

        while (start < 16) {

            // trouver prochain '>'
            int next_gt = -1;
            for (int j = start; j < 16; j++) {
                if (gt_mask & (1 << j)) {
                    next_gt = j;
                    break;
                }
            }

            int end = (next_gt == -1) ? 16 : next_gt;

            // masque segment
            uint16_t segment_mask = ((1 << end) - 1) ^ ((1 << start) - 1);
            uint16_t seg_keep = keep_mask & segment_mask;

            if (seg_keep) {
                uint8x16_t idx = vld1q_u8(shuffle_table[seg_keep]);
                uint8x16_t compacted = vqtbl1q_u8(upper, idx);

                uint8_t n = popcount_table[seg_keep];
                vst1q_u8((uint8_t*)(out + pos), compacted);
                pos += n;
                writing = 1;
            }

            // insertion newline
            if (next_gt != -1) {
                if (writing) {
                    out[pos++] = '\n';
                    writing = 0;
                }
            }

            if (next_gt == -1) break;
            start = next_gt + 1;
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
void fasta_neon_json_parser(
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
    fasta_neon_json(buffer_i, size, buffer_o, &out_size);
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

    c_time.start();
    stats_r r = stats(buffer_o, out_size);
    c_time.stop();
    std::cout << "Temps : " << c_time.ms() << " ms and debit : " << c_time.MBps(rBytes) << " Mbps" << std::endl;

    c_time.start();
    r = stats2(buffer_o, out_size);
    c_time.stop();
    std::cout << "Temps : " << c_time.ms() << " ms and debit : " << c_time.MBps(rBytes) << " Mbps" << std::endl;
    printf("A    = %10d\n", r.A);
    printf("C    = %10d\n", r.C);
    printf("T    = %10d\n", r.T);
    printf("G    = %10d\n", r.G);
    printf("size = %10zu\n", out_size);
    printf("NL   = %10d\n", r.nL);

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