#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdint>
//
#include "../tools/TimeMeasure.hpp"
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
static inline uint8_t is_valid_mask(uint8_t c) {
    c = (uint8_t)toupper(c);
    return (c=='A') | (c=='C') | (c=='G') |
           (c=='T') | (c=='U') | (c=='N');
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
void fasta_branchless(
    const char *input,
    size_t size,
    char *output,
    size_t *out_size
) {
    size_t pos = 0;

    uint8_t in_header = 0;
    uint8_t writing   = 0;

    for (size_t i = 0; i < size; i++) {
        uint8_t c = (uint8_t)input[i];

        uint8_t is_gt = (c == '>');
        uint8_t is_nl = (c == '\n');
        uint8_t is_ws = (c == ' ') | (c == '\t') | (c == '\r') | is_nl;

        // transitions header
        uint8_t end_header = in_header & is_nl;
        in_header = (in_header & ~end_header) | is_gt;

        // conditions
        uint8_t can_write = (~in_header) & (~is_ws);
        uint8_t valid = is_valid_mask(c);
        uint8_t write_mask = can_write & valid;

        // majuscule branchless
        uint8_t upper = (uint8_t)toupper(c);

        // écriture conditionnelle branchless
        output[pos] = upper;
        pos += write_mask;

        // gérer saut de ligne entre séquences
        uint8_t new_seq = is_gt & writing;

        output[pos] = '\n';
        pos += new_seq;

        // update writing (branchless)
        writing = (writing | write_mask) & (~new_seq);
    }

    // newline final
    output[pos] = '\n';
    pos += writing;

    *out_size = pos;
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
void fasta_branchless_parser(
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
    fasta_branchless(buffer_i, size, buffer_o, &out_size);
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
