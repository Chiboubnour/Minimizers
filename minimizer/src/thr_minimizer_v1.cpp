#include "functions.hpp"
#include "hls_stream.h"

//////////////////////////////////////////////////////////////////////////////////
// Thread Reader : unpack ASCII sequence from packed_sequence
//////////////////////////////////////////////////////////////////////////////////
template <int smer_size = 19>
void thread_reader(
    const uint64_t* packed_sequence,
    hls::stream<uint8_t>& stream_o
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1
#endif
    bool stop = false;
    for (int i = 0; i < 64 && !stop; i++) {
        uint64_t word = packed_sequence[i];
        for (int j = 0; j < 8; j++) {
            char c = word & 0xFF;
         
            stream_o.write(c);
            word >>= 8;
        }
    }
    // Terminate the stream
    stream_o.write(0x00);
}

template <int smer_size>
void thread_smer(
    hls::stream<uint8_t>& stream_i,
    hls::stream<ap_uint<smer_size>>& stream_o
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1
#endif
    ap_uint<smer_size> current_smer = 0;
    ap_uint<smer_size> cur_inv_smer = 0;
    constexpr int smer = smer_size / 2;
    bool first = true;

    // Build first s-mer
    for (int i = 0; i < smer - 1; i++) {
        current_smer <<= 2;
        cur_inv_smer >>= 2;
        const uint8_t nucl = stream_i.read();
        const uint8_t c_nucl = (nucl >> 1) & 0b11;
        current_smer(1, 0) = c_nucl;
        cur_inv_smer(smer_size-1, smer_size-2) = (0x2 ^ c_nucl);
    }

    while (true) {
        current_smer <<= 2;
        cur_inv_smer >>= 2;

        const uint8_t nucl = stream_i.read();
        if (nucl == 0) break;

        const uint8_t c_nucl = (nucl >> 1) & 0b11;
        current_smer(1, 0) = c_nucl;
        cur_inv_smer(smer_size-1, smer_size-2) = (0x2 ^ c_nucl);

        const ap_uint<smer_size> vmin  = min_v1<smer_size>(current_smer, cur_inv_smer);
        const ap_uint<smer_size> vhash = hash_u64<smer_size>(vmin);

        stream_o.write(vhash);
    }
    // End marker
    stream_o.write(0x00);
}




template <int window, int smer_size>
void thread_dedup(
    hls::stream<ap_uint<smer_size>>& stream_i,
    hls::stream<ap_uint<smer_size>>& stream_o
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1
    #pragma HLS ARRAY_PARTITION variable=buffer complete
#endif
    ap_uint<smer_size> buffer[window];
    ap_uint<smer_size> last_mini = -1;
    int pos = 0;
    
    // Initialize buffer with first window-1 elements
    for (int i = 0; i < window - 1; i++) {
        buffer[i] = stream_i.read();
    }
    
    // Process the rest of the stream
    while (true) {
        // Read next element
        ap_uint<smer_size> vhash = stream_i.read();
        if (vhash == 0) break;
        
        // Add new element to buffer (circular)
        buffer[(window - 1 + pos) % window] = vhash;
        
        // Find minimum in current window
        ap_uint<smer_size> min_val = buffer[pos];
        for (int i = 1; i < window; i++) {
            int idx = (pos + i) % window;
            if (buffer[idx] < min_val) {
                min_val = buffer[idx];
            }
        }
        
        // Update position
        pos = (pos + 1) % window;
        
        // Output minimizer if it's different from the last one
        if (min_val != last_mini) {
            last_mini = min_val;
            stream_o.write(min_val);
        }
    }
    
    // End marker
    stream_o.write(0x00);
}




template <int smer_size>
void thread_store(
    hls::stream<ap_uint<smer_size>>& stream_i,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nElements
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1
#endif
    int cnt = 0;
    while (true) {
        const ap_uint<smer_size> v = stream_i.read();
        if (v == 0) break;
        tab_hash[cnt++] = v;
    }
    *nElements = cnt;
}



template <int smer, int window>
void thr_minimizer(
    const uint64_t* packed_sequence,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
    constexpr int smer_size = 2 * smer;

#ifdef _PRAGMA_HLS_
    #pragma HLS DATAFLOW
    #pragma HLS INTERFACE s_axilite port=return bundle=control
    #pragma HLS INTERFACE m_axi port=packed_sequence offset=slave bundle=gmem
    #pragma HLS INTERFACE m_axi port=tab_hash offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=packed_sequence bundle=control
    #pragma HLS INTERFACE s_axilite port=tab_hash bundle=control
    #pragma HLS INTERFACE s_axilite port=nMinizrs bundle=control
#endif
    
    hls::stream<uint8_t> fifo_1("reader_to_smer");
    hls::stream<ap_uint<smer_size>> fifo_2("smer_to_dedup");
    hls::stream<ap_uint<smer_size>> fifo_3("dedup_to_store");

    thread_reader<smer_size>(packed_sequence, fifo_1);
    thread_smer<smer_size>(fifo_1, fifo_2);
    thread_dedup<window, smer_size>(fifo_2, fifo_3);
    thread_store<smer_size>(fifo_3, tab_hash, nMinizrs);
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
int t_thread_minimizer(
    const uint64_t* packed_sequence,
    const int s, const int w,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
          if(s == 28 && w == 14) { thr_minimizer<28, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 27 && w == 14) { thr_minimizer<27, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 26 && w == 14) { thr_minimizer<26, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 25 && w == 14) { thr_minimizer<25, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 24 && w == 14) { thr_minimizer<24, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 23 && w == 14) { thr_minimizer<23, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 22 && w == 14) { thr_minimizer<22, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 21 && w == 14) { thr_minimizer<21, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 20 && w == 14) { thr_minimizer<20, 14>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 19 && w == 14) { thr_minimizer<19, 14>(packed_sequence, tab_hash, nMinizrs);

    }else if(s == 28 && w == 16) { thr_minimizer<28, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 27 && w == 16) { thr_minimizer<27, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 26 && w == 16) { thr_minimizer<26, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 25 && w == 16) { thr_minimizer<25, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 24 && w == 16) { thr_minimizer<24, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 23 && w == 16) { thr_minimizer<23, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 22 && w == 16) { thr_minimizer<22, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 21 && w == 16) { thr_minimizer<21, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 20 && w == 16) { thr_minimizer<20, 16>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 19 && w == 16) { thr_minimizer<19, 16>(packed_sequence, tab_hash, nMinizrs);

    }else if(s == 28 && w == 18) { thr_minimizer<28, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 27 && w == 18) { thr_minimizer<27, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 26 && w == 18) { thr_minimizer<26, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 25 && w == 18) { thr_minimizer<25, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 24 && w == 18) { thr_minimizer<24, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 23 && w == 18) { thr_minimizer<23, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 22 && w == 18) { thr_minimizer<22, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 21 && w == 18) { thr_minimizer<21, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 20 && w == 18) { thr_minimizer<20, 18>(packed_sequence, tab_hash, nMinizrs);
    }else if(s == 19 && w == 18) { thr_minimizer<19, 18>(packed_sequence, tab_hash, nMinizrs);
        
    } else {
        if (nMinizrs) *nMinizrs = -1;
        return -1;
    }
    
    // Si nMinizrs est fourni, on y stocke le résultat
    if (nMinizrs) {
        ap_uint<64> count = 0;
        while (count < 1000000 && tab_hash[count] != 0) {
            count++;
        }
        *nMinizrs = count;
        return count;
    }
    
    // Si nMinizrs n'est pas fourni, on compte les éléments non nuls dans tab_hash
    int count = 0;
    while (count < 1000000 && tab_hash[count] != 0) {
        count++;
    }
    return count;
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
void thr_minimizer_for_synthesis(
    const uint64_t* packed_sequence,
    const int n,
    const int s,
    const int w,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
) {
    // Use the same parameters as the reference implementation
    thr_minimizer<19, 16>(packed_sequence, tab_hash, nMinizrs);
}
