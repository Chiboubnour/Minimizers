#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

#include <stdio.h>
#include <stdint.h>
#include <arm_neon.h>

#if 1
template<int W> // word length in bytes
inline uint32_t parity32(uint32_t x)
{
    x ^= x >> 16;
    x ^= x >> 8;
    x ^= x >> 4;
    x &= 0xF;
    return (0x6996 >> x) & 1;
}

template<int W> // word length in bytes
inline uint32_t parity32(uint32_t row[W])
{
    for (uint32_t x = 0; x < GF_SIZE; x += 1)
    {
        uint32_t parity = parity32(src_c & x); //
        row[x] = 1.0f - float(parity << 1);
    }
}
#endif

#define BUF 65536
//
//
//
char uppercase(const char c)
{
    return c & 0xDF;
}
//
//
//
ap_uint<8> uppercase(const ap_uint<8> c)
{
    return c & 0xDF;
}
//
//
//
template<int W> // word length in bytes
inline ap_uint<W> uppercase(const ap_uint<8 * W> value)
{
    ap_uint<W> res;
    for (int i = 0; i < W; i++)
    {
        const char lettre = value.range(8 * (i+1) - 1, 8 * i);
        res[i] = uppercase(lettre);
    }
    return vandq_u8(x,vdupq_n_u8(3));
}
//
//
//
template<int W> // word length in bytes
inline ap_uint<W> equal_char(const ap_uint<8 * W> value, const char chr)
{
    ap_uint<W> res;
    for (int i = 0; i < W; i++)
    {
        const char lettre = value.range(8 * (i+1) - 1, 8 * i);
        res[i] = (lettre == chr);
    }
    return vandq_u8(x,vdupq_n_u8(3));
}
//
//
//
template<int W> // word length in bytes
inline ap_uint<W> is_newline(const ap_uint<8 * W> value)
{
    ap_uint<W> res;
    for (int i = 0; i < W; i++)
    {
        const char lettre = value.range(8 * (i+1) - 1, 8 * i);
        res[i] = (lettre == '\n');
    }
    return vandq_u8(x,vdupq_n_u8(3));
}
//
//
//
template<int W> // word length in bytes
inline ap_uint<W> is_header_start(const ap_uint<8 * W> value)
{
    ap_uint<W> res;
    for (int i = 0; i < W; i++)
    {
        const char lettre = value.range(8 * (i+1) - 1, 8 * i);
        res[i] = (lettre == '>');
    }
    return vandq_u8(x,vdupq_n_u8(3));
}
//
//
//
/* ASCII → 2bit */
inline ap_uint<2> nucl_encode_to_2bit(const uint8_t c) {
    return (c >> 1) & 0x3;
}
//
//
//
/* packing 16 bases → 4 bytes */
template<int W> // word length in bytes
inline ap_uint<2 * W> encode_to_2bit(const ap_uint<8 * W> v)
{
    ap_uint<2 * W> res;
    for (int i = 0; i < W; i++)
    {
        const char lettre = value.range(8 * (i+1) - 1, 8 * i);
        res.range(2 * (i+1) - 1, 2 * i) = nucl_encode_to_2bit(lettre);
    }
    return res;
}
//
//
//
void parse_fastq_neon(FILE *in, FILE *out)
{
    uint8_t buf[BUF];

    int line_state = 0; /* 0..3 */
    size_t n;

    while( n = fread(buf,1,BUF,in)) )
    {
        for(size_t i=0;i+16<=n;i+=16)
        {
            uint8x16_t v =
                vld1q_u8(buf+i);

            /* newline detection */
            uint8x16_t nl =
                vceqq_u8(v,
                    vdupq_n_u8('\n'));

            if(vmaxvq_u8(nl))
                line_state =
                    (line_state+1)&3;

            /* garder uniquement ligne sequence */
            if(line_state!=1)
                continue;

            /* uppercase */
            v = vandq_u8(
                    v,
                    vdupq_n_u8(0xDF));

            /* detect bases */
            uint8x16_t isA =
                vceqq_u8(v,
                    vdupq_n_u8('A'));
            uint8x16_t isC =
                vceqq_u8(v,
                    vdupq_n_u8('C'));
            uint8x16_t isG =
                vceqq_u8(v,
                    vdupq_n_u8('G'));
            uint8x16_t isT =
                vceqq_u8(v,
                    vdupq_n_u8('T'));

            uint8x16_t valid =
                vorrq_u8(
                  vorrq_u8(isA,isC),
                  vorrq_u8(isG,isT));

            uint8x16_t code =
                ascii_to_2bit(v);

            code =
                vandq_u8(code,valid);

            uint8_t outbuf[8];
            pack16(code,outbuf);

            fwrite(outbuf,1,4,out);
        }
    }
}

/* ===============================
   ASCII → 2bit conversion
   ===============================*/
static inline uint8x16_t ascii_to_2bit(uint8x16_t v)
{
    uint8x16_t s1 = vshrq_n_u8(v,1);
    uint8x16_t s2 = vshrq_n_u8(v,2);
    uint8x16_t x  = veorq_u8(s1,s2);
    return vandq_u8(x, vdupq_n_u8(3));
}

/* ===============================
   Pack 16 bases → 4 bytes
   ===============================*/
static inline void pack16(uint8x16_t v,
                          uint8_t *out)
{
    uint8x16_t s6 = vshlq_n_u8(v,6);
    uint8x16_t s4 = vshlq_n_u8(v,4);
    uint8x16_t s2 = vshlq_n_u8(v,2);

    uint8x16_t r =
        vorrq_u8(
          vorrq_u8(s6,s4),
          vorrq_u8(s2,v));

    vst1_u8(out, vget_low_u8(r));
}

/* ===============================
   FASTA SIMD parser
   ===============================*/
template<int W> // word length in bytes
void parse_fasta(const ap_uint<8 * W>* stream_i, ap_uint<8 * W>* stream_o, ap_uint<W>* valid_o, ap_uint<W>* eos_o)
{
    uint8_t buf[BUF];

    int in_header = 0;
    size_t n;

    while((n=fread(buf,1,BUF,in)))
    {
        for(size_t i=0;i+16<=n;i+=16)
        {
            uint8x16_t v = vld1q_u8(buf+i);

            /* uppercase */
            v = vandq_u8( v, vdupq_n_u8(0xDF));

            /* detect '>' */
            uint8x16_t gt = vceqq_u8(v, vdupq_n_u8('>'));

            if(vmaxvq_u8(gt))
                in_header = 1;

            /* detect newline */
            uint8x16_t nl =
                vceqq_u8(v,
                    vdupq_n_u8('\n'));

            if(vmaxvq_u8(nl))
                in_header = 0;

            if(in_header)
                continue;

            /* detect bases */
            uint8x16_t isA   = vceqq_u8(v, vdupq_n_u8('A'));
            uint8x16_t isC   = vceqq_u8(v, vdupq_n_u8('C'));
            uint8x16_t isG   = vceqq_u8(v, vdupq_n_u8('G'));
            uint8x16_t isT   = vceqq_u8(v, vdupq_n_u8('T'));
            uint8x16_t valid = vorrq_u8( vorrq_u8(isA,isC), vorrq_u8(isG,isT));
            uint8x16_t code  = ascii_to_2bit(v);
            code = vandq_u8(code,valid);

            uint8_t outbuf[8];
            pack16(code,outbuf);
            fwrite(outbuf,1,4,out);
        }
    }
}
