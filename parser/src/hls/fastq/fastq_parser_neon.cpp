#if defined(__ARM_NEON__) || defined(__ARM_NEON)
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdint>
#include <arm_neon.h>

#include "../../tools/TimeMeasure.hpp"
#include "../../tools/stats.hpp"
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
inline uint16_t neon_movemask(const uint8x16_t input)
{
    const uint16x8_t high_bits = vreinterpretq_u16_u8 (vshrq_n_u8 (input, 7));
    const uint32x4_t paired16  = vreinterpretq_u32_u16(vsraq_n_u16(high_bits, high_bits, 7));
    const uint64x2_t paired32  = vreinterpretq_u64_u32(vsraq_n_u32(paired16, paired16, 14));
    const uint8x16_t paired64  = vreinterpretq_u8_u64 (vsraq_n_u64(paired32, paired32, 28));
    return vgetq_lane_u8(paired64, 0) | ((int)vgetq_lane_u8(paired64, 8) << 8);
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
inline int is_zero_u8x16(const uint8x16_t v) {
    return vmaxvq_u8(v) == 0;
}
//
inline int is_nzero_u8x16(const uint8x16_t v) {
    return vmaxvq_u8(v) != 0;
}
//
void fastq_neon(
    const char* input,
    size_t size,
    char* output,
    size_t* out_size
) {
    const uint8_t* ptr_i = (const uint8_t*)input;
          uint8_t* ptr_o = (      uint8_t*)output;

    const uint8x16_t n_line = vdupq_n_u8('\n');

    size_t pos = 0;
    size_t i   = 0;

    for (i = 0; i < size; i++) {

        uint16_t   m;
        uint8x16_t e;

        //
        // On saute le premier header
        //
        if( ptr_i[i] != '@' )
        {
            std::cout << "(EE) The file is not FASTQ compliant !" << std::endl;
            std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
            exit(EXIT_FAILURE);
        }

        for (; i < size; i+= 16)
        {
            const uint8x16_t v = vld1q_u8(ptr_i + i);
            e = vceqq_u8(v, n_line);
            if( is_nzero_u8x16(e) ){
                break;
            }
        }
        m  = neon_movemask(e);
        i += __builtin_ctz(m) + 1;
        
        //
        //
        //
        const int i_before = i;
        for (; i < size; i += 16)
        {
            const uint8x16_t v = vld1q_u8(ptr_i + i);
            e = vceqq_u8(v, n_line);
            vst1q_u8(ptr_o + pos, v);
            if( is_nzero_u8x16(e) )
                break;
            pos += 16;
        }
        m    = neon_movemask(e);
        i   += __builtin_ctz(m) + 1;
        pos += __builtin_ctz(m);

        const int jump_v = i - i_before - 1;

        //
        //
        //

        if( (ptr_i[i] != '+') || (ptr_i[i+1] != '\n') )
        {
            std::cout << "(EE) The file is not FASTQ compliant !" << std::endl;
            std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
            exit(EXIT_FAILURE);
        }
        i += 2;

        //
        //
        //

        i += jump_v;

        //
        //
        //

        output[pos] = '\n';
        pos += 1;
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
void fastq_neon_parser(
    const std::string& file_i,
    const std::string& file_o
){
    std::cout << "Filename: " << file_i << std::endl;

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
    fastq_neon(buffer_i, size, buffer_o, &out_size);
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