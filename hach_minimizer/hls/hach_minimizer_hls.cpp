#include <ap_int.h>
#include <hls_stream.h>

inline ap_uint<64> rotl64(ap_uint<64> x, ap_uint<8> r)
{
#pragma HLS INLINE
    return (x << r) | (x >> (64 - r));
}

inline ap_uint<64> getblock64(const ap_uint<64>* p, int i)
{
    return p[i];
}

inline ap_uint<64> fmix64(ap_uint<64> k)
{
#pragma HLS INLINE
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

inline ap_uint<64> murmur_hash_64_v1(const ap_uint<64> key, const ap_uint<64> mask)
{
#pragma HLS INLINE

    const ap_uint<64> c1 = 0x87c37b91114253d5ULL;
    const ap_uint<64> c2 = 0x4cf5ad432745937fULL;

    ap_uint<64> h1 = 42;
    ap_uint<64> h2 = 42;

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

    h1 ^= 8;
    h2 ^= 8;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    return (h1 & mask);
}

inline ap_uint<64> bfc_hash_64(ap_uint<64> key)
{
#pragma HLS INLINE
    ap_uint<64> mask = ~0ULL;

    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;

    return key;
}


void minimizer(
    const ap_uint<64>* input_minimizers,
    ap_uint<64>*       output_hashes,
    unsigned int       data_size_words
) {
#pragma HLS INTERFACE m_axi port=input_minimizers offset=slave bundle=gmem_in   max_read_burst_length=256
#pragma HLS INTERFACE m_axi port=output_hashes    offset=slave bundle=gmem_out  max_write_burst_length=256
#pragma HLS INTERFACE s_axilite port=data_size_words 
#pragma HLS INTERFACE s_axilite port=return bundle=control

main_loop:
    for (unsigned int i = 0; i < data_size_words; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT avg=100 max=100 min=100

        ap_uint<64> in  = input_minimizers[i];
        ap_uint<64> out;

#if 0
        // ===== Murmur Hash =====
        out = murmur_hash_64_v1(in, ~0ULL);
#else
        // ===== BFC Hash =====
        out = bfc_hash_64(in);
#endif

        output_hashes[i] = out;
    }
}
