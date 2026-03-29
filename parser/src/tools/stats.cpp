#include "stats.hpp"
#include <arm_neon.h>

stats_r stats(const char* memory, const size_t size)
{
    stats_r res = {0, 0, 0, 0, 0};
    for(int i = 0; i < size; i += 1) {
        const int val = memory[i];
        res.A  += (val == 'A');
        res.C  += (val == 'C');
        res.T  += (val == 'T');
        res.G  += (val == 'G');
        res.nL += (val == '\n');
    }
    return res;
}

stats_r stats2(const char* mem, const size_t size)
{
    stats_r res = {0, 0, 0, 0, 0};

    size_t i = 0;
    for (; i + 15 < size; i += 16) {
        const uint8x16_t v = vld1q_u8((const uint8_t*)(mem + i));
        const uint8x16_t A = vceqq_u8(v, vdupq_n_u8('A'));
        const uint8x16_t C = vceqq_u8(v, vdupq_n_u8('C'));
        const uint8x16_t G = vceqq_u8(v, vdupq_n_u8('G'));
        const uint8x16_t T = vceqq_u8(v, vdupq_n_u8('T'));
        const uint8x16_t n = vceqq_u8(v, vdupq_n_u8('\n'));

        res.A  += vaddvq_u8( vshrq_n_u8(A, 7) );
        res.C  += vaddvq_u8( vshrq_n_u8(C, 7) );
        res.T  += vaddvq_u8( vshrq_n_u8(T, 7) );
        res.G  += vaddvq_u8( vshrq_n_u8(G, 7) );
        res.nL += vaddvq_u8( vshrq_n_u8(n, 7) );
    }

    for (; i < size; i++) {
        const int val = mem[i];
        res.A  += (val == 'A');
        res.C  += (val == 'C');
        res.T  += (val == 'T');
        res.G  += (val == 'G');
        res.nL += (val == '\n');
    }

    return res;
}
