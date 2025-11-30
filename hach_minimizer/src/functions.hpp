#include "ap_int.h"
#include "hls_stream.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#define MEM_UNIT 64
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
inline ap_uint<2> nucl_encode_v1(const char nucl)
{
    return (nucl >> 1) & 0b11;
}
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
template<int width>
inline ap_uint<64>min_v1(const ap_uint<width> a, const ap_uint<width> b)
{
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE
#endif
    return (a < b) ? a : b;
}
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
inline ap_uint<64> mask_right(int numbits) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE
#endif
    return (numbits >= MEM_UNIT) ? ~0ULL : ((1ULL << numbits) - 1ULL);
}
//
//
//
//
// MurmurHash3 functions
//
//
//
//
#define FORCE_INLINE inline __attribute__((always_inline))
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
inline ap_uint<64> rotl64(ap_uint<64> x, ap_uint<8> r)
{
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE
#endif
    return (x << r) | (x >> (64 - r));
}
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
inline ap_uint<64> getblock64(const ap_uint<64>* p, int i)
{
    return p[i];
}
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
inline ap_uint<64> fmix64(ap_uint<64> k)
{
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE
#endif
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
template <int resu_size>
inline ap_uint<resu_size> hash_u64(uint64_t key) {
    key = (~key + (key << 21)) /*& mask*/; // key = (key << 21) - key - 1;
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8)) /*& mask*/; // key * 265
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4)) /*& mask*/; // key * 21
    key ^= key >> 28;
    key = (key + (key << 31)) /*& mask*/;
    return key;
}
//
//
//
///////////////////////////////////////////////////////////////////////////////////
//
//
//
template <int resu_size>
inline ap_uint<resu_size> murmur_hash_64_v1(const ap_uint<64> key)
{
#ifdef _PRAGMA_HLS_
    #pragma HLS PIPELINE
#endif
    const ap_uint<64> c1 = 0x87c37b91114253d5ULL;
    const ap_uint<64> c2 = 0x4cf5ad432745937fULL;

    ap_uint<64> h1 = 42; // seed blg
    ap_uint<64> h2 = 42; // seed blg

    ap_uint<64> k1 = key;
    ap_uint<64> k2 = 0;

    k1 *= c1;
    k1 = rotl64(k1, 31);
    k1 *= c2;
    h1 ^= k1;

    h1 = rotl64(h1, 27);
    h1 += h2;
    h1 = h1 * 5 + 0x52dce729;

    k2 *= c2;
    k2 = rotl64(k2, 33);
    k2 *= c1;
    h2 ^= k2;

    h2 = rotl64(h2, 31);
    h2 += h1;
    h2 = h2 * 5 + 0x38495ab5;

    h1 ^= 8; // length
    h2 ^= 8;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    return (ap_uint<resu_size>)h1;
}
//
//
//
///////////////////////////////////////////////////////////////////////////////////
//
//
//
inline ap_uint<64> murmur_hash_64_v1(const ap_uint<64> key, const ap_uint<64> mask)
{
#ifdef _PRAGMA_HLS_
    #pragma HLS PIPELINE
#endif
    const ap_uint<64> c1 = 0x87c37b91114253d5ULL;
    const ap_uint<64> c2 = 0x4cf5ad432745937fULL;

    ap_uint<64> h1 = 42; // seed blg
    ap_uint<64> h2 = 42; // seed blg

    ap_uint<64> k1 = key;
    ap_uint<64> k2 = 0;

    k1 *= c1;
    k1 = rotl64(k1, 31);
    k1 *= c2;
    h1 ^= k1;

    h1 = rotl64(h1, 27);
    h1 += h2;
    h1 = h1 * 5 + 0x52dce729;

    k2 *= c2;
    k2 = rotl64(k2, 33);
    k2 *= c1;
    h2 ^= k2;

    h2 = rotl64(h2, 31);
    h2 += h1;
    h2 = h2 * 5 + 0x38495ab5;

    h1 ^= 8; // length
    h2 ^= 8;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    return (h1 & mask);
}
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
inline void show_x_mer(char* sequence, uint64_t hash, int n, uint64_t vmin) {
    printf(" smer : ");
    for (int i = 0; i < n; i++) {
        printf("%c", sequence[i]);
    }
    printf("\n");
}
//
//
//
///////////////////////////////////////////////////////////////////////////////////
//
//
//
template<int smer_size>
inline void show_x_mer(const ap_uint<smer_size>* buffer, int n) {
    printf("buffer :");
    for (int i = 0; i < n; i++)
    {
        if( i%4 == 0 )
            printf("\n    ");
        printf("0x%10.10llX  ", buffer[i].to_uint64());
    }
    printf("\n");
}
//
//
//
///////////////////////////////////////////////////////////////////////////////////
//
//
//
