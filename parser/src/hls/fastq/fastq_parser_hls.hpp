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
extern void fastq_hls(
    const char* input,
    size_t      size,
    char*       output,
    size_t*     out_size
);
//
//
extern void fastq_neon_hls(
    const std::string& file_i,
    const std::string& file_o
);
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//

//uint8x16_t 

template<int P, int W>
struct ap_uint8xP_t
{
    ap_uint<W> v[P];
    ap_uint<W> valid;
    ap_uint<W> eos;
    ap_uint<W> eof];
};