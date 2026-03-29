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
void fastq_branchless(
    const char *input,
    size_t size,
    char *output,
    size_t *out_size
) {
    size_t pos = 0;
    size_t i = 0;
    for (i = 0; i < size; i++) {

        if( input[i] != '@' )
        {
            printf("Error @ (char = %c, pos = %zu)\n", input[i], pos);
            exit( EXIT_FAILURE );
        }

        //
        // On saute le premier header
        //
        while (input[i] != '\n')
        {
            i += 1;
        }
        i += 1; // on skip maintenant le '\n'
        
        while (input[i] != '\n')
        {
            output[pos] = input[i];
            i   += 1;
            pos += 1;
        }
        i += 1; // on skip maintenant le '\n'
        
        if( input[i] != '+' )
        {
            printf("Error + (char = %c, pos = %zu)\n", input[i], pos);
            exit( EXIT_FAILURE );
        }

        while (input[i] != '\n')
        {
            printf("%c", input[i]); fflush( stdout );
            i += 1;
        }
        i += 1; // on skip maintenant le '\n'

        while( (input[i] != '\n') && (i < size) )
        {
            i += 1;
        }
        // i += 1; non car deja rajouté dans la boucle FOR
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
void fastq_branchless_parser(
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
    fastq_branchless(buffer_i, size, buffer_o, &out_size);
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
