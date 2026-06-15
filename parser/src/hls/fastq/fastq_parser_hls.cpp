#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdint>

#include "fastq_parser_hls.hpp"
#include "../../tools/TimeMeasure.hpp"
#include "../../tools/stats.hpp"
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
bool isEqual(const ap_uint<8> val, const uint8_t c)
{
    return (val == c);
}
//
template<int P>
bool has_new_line(const ap_uint<8*P>& lane)
{
    bool hasNL = false;
    for(int i = 0; i < P; i += 1){
        const ap_uint<8> val = lane.range(8 * i - 1, 8 * i);
        hasNL |= isEqual(v, '\n');
    }
    return hasNL;
}
//
template<int P, int W>
void fastq_hls(
    const ap_uint<P>* input,
    size_t size,
    ap_uintWxP_t<P,W> output
) {
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
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
void fastq_hls_parser(
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
    fastq_hls(buffer_i, size, buffer_o);
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
