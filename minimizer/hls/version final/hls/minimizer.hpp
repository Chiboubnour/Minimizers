#ifndef MINIMIZER_64BITS_2_PARAMS_HPP
#define MINIMIZER_64BITS_2_PARAMS_HPP

constexpr int PAR_FACTOR  = 8;
constexpr int SMER_SIZE   = 21;
constexpr int WINDOW_SIZE = 11;
constexpr int SMER_BITS   = 2 * SMER_SIZE;
constexpr int MEM_WIDTH   = 64;

static_assert(MEM_WIDTH % (8 * PAR_FACTOR) == 0,
              "MEM_WIDTH doit etre un multiple de 8*PAR_FACTOR");
constexpr int MEM_RATIO = MEM_WIDTH / (8 * PAR_FACTOR);


constexpr int MAX_BASES_HW = 1000000;


constexpr int COSIM_MAX_BASES = MAX_BASES_HW;
constexpr int COSIM_MAX_WORDS = (COSIM_MAX_BASES + PAR_FACTOR - 1) / PAR_FACTOR;
constexpr int COSIM_MAX_BEATS = (COSIM_MAX_WORDS + MEM_RATIO - 1) / MEM_RATIO;

constexpr int GROUP_W       = 2;
constexpr int GROUP_BITS    = 64 * GROUP_W;
constexpr int COSIM_MAX_HASHGROUP = (COSIM_MAX_BASES + GROUP_W - 1) / GROUP_W;

#endif
