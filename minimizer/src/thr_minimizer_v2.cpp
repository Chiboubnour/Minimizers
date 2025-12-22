#include "functions.hpp"
#include <hls_stream.h>
#include <ap_int.h>

//
//
//
///////////////////////////////////////////////////////////////////////////////////
//
//
//
////#define _DEBUG_
//
//
template <int smer_size = 19>
void thread_reader_v2(
    const uint64_t* packed_sequence,
    hls::stream<uint8_t>& stream_o
) {
    for (int i = 0; i < 64; i++) {
        uint64_t word = packed_sequence[i];
        for (int j = 0; j < 8; j++) {
            uint8_t c = word & 0xFF;
            if (c == 0) { stream_o.write(0x00); return; }
            stream_o.write(c);
            word >>= 8;
        }
    }
    stream_o.write(0x00); // fin de séquence
}

template <int smer_size>
void thread_smer_v2(
    hls::stream<uint8_t>& stream_i,
    hls::stream<ap_uint<smer_size>>& stream_o
) {
    ap_uint<smer_size> current_smer = 0;
    ap_uint<smer_size> cur_inv_smer = 0;
    constexpr int smer = smer_size / 2;

    // Build first s-mer
    for (int i = 0; i < smer - 1; i++) {
        current_smer <<= 2;
        cur_inv_smer >>= 2;
        uint8_t nucl = stream_i.read();
        if (nucl == 0) { stream_o.write(0x00); return; }
        uint8_t c_nucl = (nucl >> 1) & 0x03;
        current_smer(1,0) = c_nucl;
        cur_inv_smer(smer_size-1, smer_size-2) = 0x2 ^ c_nucl;
    }

    while (true) {
        current_smer <<= 2;
        cur_inv_smer >>= 2;
        uint8_t nucl = stream_i.read();
        if (nucl == 0) break;
        uint8_t c_nucl = (nucl >> 1) & 0x03;
        current_smer(1,0) = c_nucl;
        cur_inv_smer(smer_size-1, smer_size-2) = 0x2 ^ c_nucl;

        ap_uint<smer_size> vmin  = min_v1<smer_size>(current_smer, cur_inv_smer);
        ap_uint<smer_size> vhash = hash_u64<smer_size>(vmin);
        stream_o.write(vhash);
    }
    stream_o.write(0x00);
}

template <int window, int smer_size>
void thread_dedup_v2(
    hls::stream<ap_uint<smer_size>>& stream_i,
    hls::stream<ap_uint<smer_size>>& stream_o
) {
    ap_uint<smer_size> buffer[window];
    ap_uint<smer_size> last_mini = -1;

    for (int i = 0; i < window-1; i++) buffer[i] = stream_i.read();

    int pos = 0;
    while (true) {
        ap_uint<smer_size> vhash = stream_i.read();
        if (vhash == 0) break;

        buffer[(window-1+pos)%window] = vhash;

        ap_uint<smer_size> min_val = buffer[pos];
        for (int i = 1; i < window; i++) {
            int idx = (pos+i)%window;
            if (buffer[idx] < min_val) min_val = buffer[idx];
        }

        pos = (pos+1)%window;

        if (min_val != last_mini) {
            last_mini = min_val;
            stream_o.write(min_val);
        }
    }
    stream_o.write(0x00);
}

template <int smer_size>
void thread_store_v2(
    hls::stream<ap_uint<smer_size>>& stream_i,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nElements
) {
    int cnt = 0;
    while (true) {
        ap_uint<smer_size> vHash = stream_i.read();
        if (vHash == 0) break;
        tab_hash[cnt++] = vHash;
    }
    *nElements = cnt;
}

//
//
//
///////////////////////////////////////////////////////////////////////////////////
//
//
//
#define DATA_DEPTH 64
//
//
//
template <int smer, int window>
void thr_minimizer_v2(
    const uint64_t* packed_sequence,
    ap_uint<64> n,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INTERFACE mode=m_axi     port=packed_sequence
    #pragma HLS INTERFACE mode=s_axilite port=n
    #pragma HLS INTERFACE mode=m_axi     port=nMinizrs
    #pragma HLS INTERFACE mode=s_axilite port=return bundle=control
    #pragma HLS DATAFLOW
#endif
    constexpr int smer_size = 2 * smer;

    hls::stream< uint8_t, DATA_DEPTH > fifo_1;
    hls::stream< ap_uint< smer_size >, DATA_DEPTH > fifo_2;
    hls::stream< ap_uint< smer_size >, DATA_DEPTH > fifo_3;
    hls::stream< ap_uint< smer_size >, DATA_DEPTH > fifo_4;

#ifdef _DEBUG_
    char* kmer_str = (char*)packed_sequence;
    printf("Sequence (1) = %s\n", kmer_str);
#endif

thread_reader_v2<smer_size>(packed_sequence, fifo_1);
thread_smer_v2<smer_size>(fifo_1, fifo_2);
thread_dedup_v2<window, smer_size>(fifo_2, fifo_3);
thread_store_v2<smer_size>(fifo_3, tab_hash, nMinizrs);


#ifdef _DEBUG_
    printf("#elements in fifo_1: %llu\n", fifo_1.size());
    printf("#elements in fifo_2: %llu\n", fifo_2.size());
    printf("#elements in fifo_3: %llu\n", fifo_3.size());
#endif
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
void t_thread_minimizer_v2(
    const uint64_t* packed_sequence,
    const int s, const int w,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
          if(s == 28 && w == 14) { thr_minimizer_v2<28, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 27 && w == 14) { thr_minimizer_v2<27, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 26 && w == 14) { thr_minimizer_v2<26, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 25 && w == 14) { thr_minimizer_v2<25, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 24 && w == 14) { thr_minimizer_v2<24, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 23 && w == 14) { thr_minimizer_v2<23, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 22 && w == 14) { thr_minimizer_v2<22, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 21 && w == 14) { thr_minimizer_v2<21, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 20 && w == 14) { thr_minimizer_v2<20, 14>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 19 && w == 14) { thr_minimizer_v2<19, 14>(packed_sequence, s, tab_hash, nMinizrs);

    }else if(s == 28 && w == 16) { thr_minimizer_v2<28, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 27 && w == 16) { thr_minimizer_v2<27, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 26 && w == 16) { thr_minimizer_v2<26, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 25 && w == 16) { thr_minimizer_v2<25, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 24 && w == 16) { thr_minimizer_v2<24, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 23 && w == 16) { thr_minimizer_v2<23, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 22 && w == 16) { thr_minimizer_v2<22, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 21 && w == 16) { thr_minimizer_v2<21, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 20 && w == 16) { thr_minimizer_v2<20, 16>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 19 && w == 16) { thr_minimizer_v2<19, 16>(packed_sequence, s, tab_hash, nMinizrs);

    }else if(s == 28 && w == 18) { thr_minimizer_v2<28, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 27 && w == 18) { thr_minimizer_v2<27, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 26 && w == 18) { thr_minimizer_v2<26, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 25 && w == 18) { thr_minimizer_v2<25, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 24 && w == 18) { thr_minimizer_v2<24, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 23 && w == 18) { thr_minimizer_v2<23, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 22 && w == 18) { thr_minimizer_v2<22, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 21 && w == 18) { thr_minimizer_v2<21, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 20 && w == 18) { thr_minimizer_v2<20, 18>(packed_sequence, s, tab_hash, nMinizrs);
    }else if(s == 19 && w == 18) { thr_minimizer_v2<19, 18>(packed_sequence, s, tab_hash, nMinizrs);
    } else {
#ifndef _PRAGMA_HLS_
        std::cerr << "Error: Unsupported s value: " << s << std::endl;
        std::cerr << "Supported values are between 16 and 31." << std::endl;
        std::cerr << "Please check the input parameters." << std::endl;
        std::cerr << "Exiting..." << std::endl;
#endif
        *nMinizrs = -1;
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
void thr_minimizer_for_synthesis_v2(
    const uint64_t* packed_sequence,
    const int s,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
        thr_minimizer_v2<19, 16>(packed_sequence, s, tab_hash, nMinizrs);
}
