#if defined(__ARM_NEON__) || defined(__ARM_NEON)
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include <cstdint>
#include <string>
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
extern void fasta_neon_json(
    const char* input,
    size_t      size,
    char*       output,
    size_t*     out_size
);
//
//
extern void fasta_neon_json_parser(
    const std::string& file_i,
    const std::string& file_o
);
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
#endif