#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdint>
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
extern void fastq_neon(
    const char* input,
    size_t      size,
    char*       output,
    size_t*     out_size
);
//
//
extern void fastq_neon_parser(
    const std::string& file_i,
    const std::string& file_o
);
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
