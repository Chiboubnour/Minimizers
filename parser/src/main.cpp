#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cinttypes>
#include <ap_int.h>
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
#include "./hls/fastq/fastq_parser.hpp"
#include "./hls/fastq/fastq_parser_neon.hpp"
//
#include "./hls/fasta/fasta_parser.hpp"
#include "./hls/fasta/fasta_parser_neon.hpp"
#include "./hls/fasta/fasta_parser_neon_ultimate.hpp"
#include "./hls/fasta/fasta_parser_neon_json.hpp"
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
int main(int argc, char* argv[])
{
    if( argc < 4 ){
        printf("(II) Not enought arguments\n");
        printf("(II) ./main input_file output_file type\n");
        exit(EXIT_FAILURE);
    }

    std::string file_i = argv[1];
    std::string file_o = argv[2];
    std::string type   = argv[3];

    if( type == "naive" ){
        fasta_branchless_parser(file_i, file_o);
#if defined(__ARM_NEON__) || defined(__ARM_NEON)
    }else if( type == "neon" ){
        fasta_neon_parser(file_i, file_o);
    }else if( type == "ultra" ){
        fasta_neon_ultimate_parser(file_i, file_o);
    }else if( type == "json" ){
        fasta_neon_json_parser(file_i, file_o);
#endif
    }else if( type == "naive_q" ){
        fastq_branchless_parser(file_i, file_o);
#if defined(__ARM_NEON__) || defined(__ARM_NEON)
    }else if( type == "neon_q" ){
        fastq_neon_parser(file_i, file_o);
#endif
    }else{
        printf("(II) Parser type not recognize (%s)\n", type.c_str());
        exit(EXIT_FAILURE);
    }
    
    return 0;
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
