#include "functions.hpp"
//
//
//
///////////////////////////////////////////////////////////////////////////////////
//
//
//
template <int smer, int window>
int t_minimizer_sequence(
    const uint64_t* packed_sequence,
    ap_uint<64> n,
    ap_uint<64>* tab_hash//,
//    ap_uint<64>* tab_occu
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INTERFACE mode=m_axi     port=packed_sequence
    #pragma HLS INTERFACE mode=s_axilite port=n
    #pragma HLS INTERFACE mode=m_axi     port=tab_hash
    #pragma HLS INTERFACE mode=s_axilite port=return      bundle=control
#endif
    constexpr int smer_size = 2 * smer;
#if 0
    printf("=> sequence %s\n", (char*)packed_sequence);
#endif
    //
    // On realise la formation du premier kmer
    //

    ap_uint<smer_size> current_smer = 0;
    ap_uint<smer_size> cur_inv_smer = 0;

    int cnt = 0;
    //
    //
    ////////////////////////////////////////////////////////////////////////////////////
    //
    //
    bool stop = false;
    int n_elements = 0;
    char sequence[512];
    for (int i = 0; i < 64; i++)
    {
#ifdef _PRAGMA_HLS_
    #pragma HLS PIPELINE
#endif
        ap_uint<64> word = packed_sequence[i];
        for (int j = 0; j < 8; ++j)
        {
            if( word.to_uchar() == 0 ){ // caractere nul indique la fin de la sequence
                stop = true;            // de nucleotides
                break;
            }
            sequence[n_elements++] = word;
            word                   = word >> 8;
        }
        if( stop ) break;
    }
#if 0
    printf("=> sequence size %d\n", n_elements);
#endif 

    //
    //
    ////////////////////////////////////////////////////////////////////////////////////
    //
    //
    ap_uint<smer_size> buffer[window]; // nombre de m-mer/k-mer
    for(int x = 0; x < window; x += 1)
        buffer[x] = 0;


    //
    //
    ////////////////////////////////////////////////////////////////////////////////////
    //
    //
    for (int i = 0; i < smer - 1; i++)
    {
#ifdef _PRAGMA_HLS_
    #pragma HLS PIPELINE
#endif
        current_smer <<= 2;
        cur_inv_smer >>= 2;
        const ap_uint<8>   nucl = sequence[i];                   // extract the ASCII char
        const ap_uint<2> c_nucl = (nucl >> 1) & 0b11;            // convert to 2 bits representation
        current_smer(          1,           0) = c_nucl;         // update the normal  smer
        cur_inv_smer(smer_size-1, smer_size-2) = (0x2 ^ c_nucl); // update the inverse smer
    }

    ap_uint<smer_size>  last_mini = -1;
    for(int i = smer - 1; i < n_elements; i += 1)
    {
#ifdef _PRAGMA_HLS_
    #pragma HLS PIPELINE
#endif
        current_smer <<= 2;
        cur_inv_smer >>= 2;
        const ap_uint<8>   nucl = sequence[i];                   // extract the ASCII char
        const ap_uint<2> c_nucl = (nucl >> 1) & 0b11;            // convert to 2 bits representation
        current_smer(          1,           0) = c_nucl;         // update the normal  smer
        cur_inv_smer(smer_size-1, smer_size-2) = (0x2 ^ c_nucl); // update the inverse smer

        const ap_uint<smer_size> vmin  = min_v1<smer_size>( current_smer, cur_inv_smer);
        const ap_uint<smer_size> vhash = hash_u64/*murmur_hash_64_v1*/<smer_size>( vmin );
//      printf("-> min(%8.8llX, %8.8llX) = %8.8llX  => hash = 0x%16.16llX\n", current_smer.to_uint64(), cur_inv_smer.to_uint64(), vmin.to_uint64(), vhash.to_uint64());
#if 0
        show_x_mer(sequence + i - smer + 1, vhash.to_uint64(), smer, vmin.to_uint64());
        printf("  - min(%8.8llX, %8.8llX) = %8.8llX  => hash = 0x%16.16llX\n", current_smer.to_uint64(), cur_inv_smer.to_uint64(), vmin.to_uint64(), vhash.to_uint64());
#endif

        //
        // On calcule le minimiseur de la fenetre Z
        //
#if 0        
        show_x_mer<smer_size>(buffer, z);
#endif
        ap_uint<64> minz = vhash;
        for(int p = 0; p < window; p += 1) {
            if( buffer[p] < minz )
                minz = buffer[p];
        }
#if 0
        show_x_mer(sequence + i - smer + 1, minz.to_uint64(), smer, vmin.to_uint64());
#endif
        //
        // On fait vieillir la fenetre des valeurs de hash
        //

        for(int p = 0; p < window - 1; p += 1) {
            buffer[p] = buffer[p+1];
        }
        buffer[window-1] = vhash;

#if 0
        show_x_mer<smer_size>(buffer, z);
#endif
        //
        // On memorise apres avoir dé-dupliqué
        //

        if( ((i-smer+1) >= window) && (minz != last_mini) )

        {
#if 0
            printf("=> storing s-mer [0x%16.16llX] because last_mini = [0x%16.16llX]\n", minz.to_uint64(), last_mini.to_uint64());
#endif
            last_mini        = minz;
            tab_hash [cnt]   = minz;
//          tab_occu [cnt]   = 1;
            cnt++;
#if 0
            bool a = ((i-smer+1) >= z);
            bool b = (minz != last_mini);
            printf("a = %d et b = %d\n", a, b);
#endif
        }
        else
        {
#if 0
            printf("=> SKIPPING s-mer [0x%16.16llX] because last_mini = [0x%16.16llX]\n", minz.to_uint64(), last_mini.to_uint64());
#endif
#if 0
            bool a = ((i-smer+1) >= z);
            bool b = (minz != last_mini);
            printf("a = %d et b = %d\n", a, b);
#endif
//          if(cnt != 0)
//              tab_occu[cnt-1] += 1;
        }
    }

    return cnt;
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
int t_minimizer_sequence(
    const uint64_t* packed_sequence,
    const int s,
    const int w,
    ap_uint<64>* tab_hash
) {
          if(s == 28 && w == 14) { return t_minimizer_sequence<28, 14>(packed_sequence, s, tab_hash);
    }else if(s == 27 && w == 14) { return t_minimizer_sequence<27, 14>(packed_sequence, s, tab_hash);
    }else if(s == 26 && w == 14) { return t_minimizer_sequence<26, 14>(packed_sequence, s, tab_hash);
    }else if(s == 25 && w == 14) { return t_minimizer_sequence<25, 14>(packed_sequence, s, tab_hash);
    }else if(s == 24 && w == 14) { return t_minimizer_sequence<24, 14>(packed_sequence, s, tab_hash);
    }else if(s == 23 && w == 14) { return t_minimizer_sequence<23, 14>(packed_sequence, s, tab_hash);
    }else if(s == 22 && w == 14) { return t_minimizer_sequence<22, 14>(packed_sequence, s, tab_hash);
    }else if(s == 21 && w == 14) { return t_minimizer_sequence<21, 14>(packed_sequence, s, tab_hash);
    }else if(s == 20 && w == 14) { return t_minimizer_sequence<20, 14>(packed_sequence, s, tab_hash);
    }else if(s == 19 && w == 14) { return t_minimizer_sequence<19, 14>(packed_sequence, s, tab_hash);
    }else if(s == 28 && w == 16) { return t_minimizer_sequence<28, 16>(packed_sequence, s, tab_hash);
    }else if(s == 27 && w == 16) { return t_minimizer_sequence<27, 16>(packed_sequence, s, tab_hash);
    }else if(s == 26 && w == 16) { return t_minimizer_sequence<26, 16>(packed_sequence, s, tab_hash);
    }else if(s == 25 && w == 16) { return t_minimizer_sequence<25, 16>(packed_sequence, s, tab_hash);
    }else if(s == 24 && w == 16) { return t_minimizer_sequence<24, 16>(packed_sequence, s, tab_hash);
    }else if(s == 23 && w == 16) { return t_minimizer_sequence<23, 16>(packed_sequence, s, tab_hash);
    }else if(s == 22 && w == 16) { return t_minimizer_sequence<22, 16>(packed_sequence, s, tab_hash);
    }else if(s == 21 && w == 16) { return t_minimizer_sequence<21, 16>(packed_sequence, s, tab_hash);
    }else if(s == 20 && w == 16) { return t_minimizer_sequence<20, 16>(packed_sequence, s, tab_hash);
    }else if(s == 19 && w == 16) { return t_minimizer_sequence<19, 16>(packed_sequence, s, tab_hash);
    }else if(s == 28 && w == 18) { return t_minimizer_sequence<28, 18>(packed_sequence, s, tab_hash);
    }else if(s == 27 && w == 18) { return t_minimizer_sequence<27, 18>(packed_sequence, s, tab_hash);
    }else if(s == 26 && w == 18) { return t_minimizer_sequence<26, 18>(packed_sequence, s, tab_hash);
    }else if(s == 25 && w == 18) { return t_minimizer_sequence<25, 18>(packed_sequence, s, tab_hash);
    }else if(s == 24 && w == 18) { return t_minimizer_sequence<24, 18>(packed_sequence, s, tab_hash);
    }else if(s == 23 && w == 18) { return t_minimizer_sequence<23, 18>(packed_sequence, s, tab_hash);
    }else if(s == 22 && w == 18) { return t_minimizer_sequence<22, 18>(packed_sequence, s, tab_hash);
    }else if(s == 21 && w == 18) { return t_minimizer_sequence<21, 18>(packed_sequence, s, tab_hash);
    }else if(s == 20 && w == 18) { return t_minimizer_sequence<20, 18>(packed_sequence, s, tab_hash);
    }else if(s == 19 && w == 18) { return t_minimizer_sequence<19, 18>(packed_sequence, s, tab_hash);

    } else {
//#ifdef _DEBUG_
        std::cerr << "Error: Unsupported s value: " << s << std::endl;
        std::cerr << "Error: Unsupported s value: " << w << std::endl;
        std::cerr << "Supported values are between 16 and 31." << std::endl;
        std::cerr << "Please check the input parameters." << std::endl;
        std::cerr << "Exiting..." << std::endl;
//#endif
        return -1;
    }
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
int t_minimizer_v1_for_synthesis(
    const uint64_t* packed_sequence,
    const int s,
    ap_uint<64>* tab_hash
) {
        return t_minimizer_sequence<7,4>(packed_sequence, s, tab_hash/*, tab_occu*/);
}