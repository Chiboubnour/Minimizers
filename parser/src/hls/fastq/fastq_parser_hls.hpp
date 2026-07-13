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
template<int P, int W>
struct ap_uintWxP_t
{
    ap_uint<W> v[P];    // the ASCII values
    ap_uint<P> valid;   // the valid bit : is it sequence information
    ap_uint<P> eos;     // end of sequence bit
    ap_uint<P> eof;     // end of file bit => we should stop
};

template <int gf_size>
struct symbols_i
{
	int32_t value[gf_size];
	bool is_freq;
};
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
template<int P, int W>
extern void fastq_hls(
    const char* input,
    size_t      size,
    ap_uintWxP_t<P,W>* output
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
