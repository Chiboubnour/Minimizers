#include "functions.hpp"

//////////////////////////////////////////////////////////////////////////////////////
/// VERSION TEMPLATE (VRAIE IMPLÉMENTATION)
//////////////////////////////////////////////////////////////////////////////////////

template <int smer, int window>
int t_minimizer_v1(
    const uint64_t* packed_sequence,
    const int n,
    ap_uint<64>* tab_hash
) {
    constexpr int smer_size = 2 * smer;

    ap_uint<smer_size> current_smer = 0;
    ap_uint<smer_size> cur_inv_smer = 0;

    int cnt = 0;

    //////////////////////////////////////////////////////////////////////////
    // 1. Dépaquetage : packed_sequence -> sequence ASCII
    //////////////////////////////////////////////////////////////////////////

    bool stop = false;
    int seq_len = 0;
    char sequence[512];

    for (int i = 0; i < 64; i++)
    {
        ap_uint<64> word = packed_sequence[i];

        for (int j = 0; j < 8; ++j)
        {
            if(seq_len >= 512){
                stop = true;
                break;
            }

            char c = word & 0xFF; // octet de droite

            if(c == 0){
                stop = true;
                break;
            }

            sequence[seq_len++] = c;
            word >>= 8;
        }

        if(stop) break;
    }

    //////////////////////////////////////////////////////////////////////////
    // 2. Initialisation buffer de fenêtre
    //////////////////////////////////////////////////////////////////////////

    ap_uint<smer_size> buffer[window];
    for(int x = 0; x < window; x++) {
        buffer[x] = 0;
    }

    //////////////////////////////////////////////////////////////////////////
    // 3. Construction du premier s-mer
    //////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < smer - 1; i++)
    {
        current_smer <<= 2;
        cur_inv_smer >>= 2;

        const ap_uint<8>   nucl   = sequence[i];
        const ap_uint<2>   c_nucl = (nucl >> 1) & 0b11;

        current_smer(1,0) = c_nucl;
        cur_inv_smer(smer_size-1, smer_size-2) = (0x2 ^ c_nucl);
    }

    //////////////////////////////////////////////////////////////////////////
    // 4. Parcours + calcul minimizers
    //////////////////////////////////////////////////////////////////////////

    ap_uint<smer_size> last_mini = -1;

    for(int i = smer - 1; i < seq_len; i++)
    {
        current_smer <<= 2;
        cur_inv_smer >>= 2;

        const ap_uint<8>   nucl   = sequence[i];
        const ap_uint<2>   c_nucl = (nucl >> 1) & 0b11;

        current_smer(1,0) = c_nucl;
        cur_inv_smer(smer_size-1, smer_size-2) = (0x2 ^ c_nucl);

        const ap_uint<smer_size> vmin  = min_v1<smer_size>(current_smer, cur_inv_smer);
        const ap_uint<smer_size> vhash = hash_u64<smer_size>(vmin);

        ap_uint<64> minz = vhash;

        for(int p = 0; p < window; p++) {
            if(buffer[p] < minz)
                minz = buffer[p];
        }

        for(int p = 0; p < window - 1; p++) {
            buffer[p] = buffer[p+1];
        }

        buffer[window - 1] = vhash;

        if(((i - smer + 1) >= window) && (minz != last_mini))
        {
            last_mini = minz;
            tab_hash[cnt++] = minz;
        }
    }

    return cnt;
}

//////////////////////////////////////////////////////////////////////////////////////
/// VERSION GENERIC — Dispatcher s/w
//////////////////////////////////////////////////////////////////////////////////////

int t_minimizer_sequence(
    const uint64_t* packed_sequence,
    const int n,
    const int s,
    const int w,
    ap_uint<64>* tab_hash
) {
         if(s == 28 && w == 14) return t_minimizer_v1<28,14>(packed_sequence, n, tab_hash);
    else if(s == 27 && w == 14) return t_minimizer_v1<27,14>(packed_sequence, n, tab_hash);
    else if(s == 26 && w == 14) return t_minimizer_v1<26,14>(packed_sequence, n, tab_hash);
    else if(s == 25 && w == 14) return t_minimizer_v1<25,14>(packed_sequence, n, tab_hash);
    else if(s == 24 && w == 14) return t_minimizer_v1<24,14>(packed_sequence, n, tab_hash);
    else if(s == 23 && w == 14) return t_minimizer_v1<23,14>(packed_sequence, n, tab_hash);
    else if(s == 22 && w == 14) return t_minimizer_v1<22,14>(packed_sequence, n, tab_hash);
    else if(s == 21 && w == 14) return t_minimizer_v1<21,14>(packed_sequence, n, tab_hash);
    else if(s == 20 && w == 14) return t_minimizer_v1<20,14>(packed_sequence, n, tab_hash);
    else if(s == 19 && w == 14) return t_minimizer_v1<19,14>(packed_sequence, n, tab_hash);

    else if(s == 28 && w == 16) return t_minimizer_v1<28,16>(packed_sequence, n, tab_hash);
    else if(s == 27 && w == 16) return t_minimizer_v1<27,16>(packed_sequence, n, tab_hash);
    else if(s == 26 && w == 16) return t_minimizer_v1<26,16>(packed_sequence, n, tab_hash);
    else if(s == 25 && w == 16) return t_minimizer_v1<25,16>(packed_sequence, n, tab_hash);
    else if(s == 24 && w == 16) return t_minimizer_v1<24,16>(packed_sequence, n, tab_hash);
    else if(s == 23 && w == 16) return t_minimizer_v1<23,16>(packed_sequence, n, tab_hash);
    else if(s == 22 && w == 16) return t_minimizer_v1<22,16>(packed_sequence, n, tab_hash);
    else if(s == 21 && w == 16) return t_minimizer_v1<21,16>(packed_sequence, n, tab_hash);
    else if(s == 20 && w == 16) return t_minimizer_v1<20,16>(packed_sequence, n, tab_hash);
    else if(s == 19 && w == 16) return t_minimizer_v1<19,16>(packed_sequence, n, tab_hash);

    else if(s == 28 && w == 18) return t_minimizer_v1<28,18>(packed_sequence, n, tab_hash);
    else if(s == 27 && w == 18) return t_minimizer_v1<27,18>(packed_sequence, n, tab_hash);
    else if(s == 26 && w == 18) return t_minimizer_v1<26,18>(packed_sequence, n, tab_hash);
    else if(s == 25 && w == 18) return t_minimizer_v1<25,18>(packed_sequence, n, tab_hash);
    else if(s == 24 && w == 18) return t_minimizer_v1<24,18>(packed_sequence, n, tab_hash);
    else if(s == 23 && w == 18) return t_minimizer_v1<23,18>(packed_sequence, n, tab_hash);
    else if(s == 22 && w == 18) return t_minimizer_v1<22,18>(packed_sequence, n, tab_hash);
    else if(s == 21 && w == 18) return t_minimizer_v1<21,18>(packed_sequence, n, tab_hash);
    else if(s == 20 && w == 18) return t_minimizer_v1<20,18>(packed_sequence, n, tab_hash);
    else if(s == 19 && w == 18) return t_minimizer_v1<19,18>(packed_sequence, n, tab_hash);

    else return -1;
}

//////////////////////////////////////////////////////////////////////////////////////
/// VERSION POUR LA SYNTHESE
//////////////////////////////////////////////////////////////////////////////////////

int t_minimizer_v1_for_synthesis(
    const uint64_t* packed_sequence,
    const int n,
    const int s,
    const int w,
    ap_uint<64>* tab_hash
) {
    return t_minimizer_v1<31,19>(packed_sequence, n, tab_hash);
}
