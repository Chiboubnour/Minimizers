#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

#define SMER        19
#define SMER_SIZE   (2 * SMER)
#define WINDOW_SIZE 16
#define DATA_DEPTH  1024

inline ap_uint<2> nucl_encode(uint8_t c) {
#pragma HLS INLINE
    return (c >> 1) & 0x3;
}


inline ap_uint<SMER_SIZE> min_v1(ap_uint<SMER_SIZE> a, ap_uint<SMER_SIZE> b) {
#pragma HLS INLINE
    return (a < b) ? a : b;
}


inline ap_uint<64> hash_u64(ap_uint<SMER_SIZE> key) {
    ap_uint<64> k64 = key;
    k64 = ~k64 + (k64 << 21);
    k64 ^= k64 >> 24;
    k64 = ((k64 + (k64 << 3)) + (k64 << 8));
    k64 ^= k64 >> 14;
    k64 = ((k64 + (k64 << 2)) + (k64 << 4));
    k64 ^= k64 >> 28;
    k64 = k64 + (k64 << 31);
    return k64;
}

void thread_reader(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n,
    hls::stream<uint8_t>& stream_o
) {
    ap_uint<64> nb_words = (n + 7) >> 3;
    for (ap_uint<64> i = 0; i < nb_words; i++) {
#pragma HLS PIPELINE II=1
        ap_uint<64> word = packed_sequence[i];

        for (int j = 0; j < 8; j++) {
#pragma HLS UNROLL
            uint64_t base_idx = (i << 3) + j;
            if (base_idx >= n) break;
            uint8_t c = word.range(8*j + 7, 8*j);
            stream_o.write(c);
        }
    }
    stream_o.write(0);
}


// ============================================================
// THREAD S-MER
// ============================================================
void thread_smer(
    hls::stream<uint8_t>& stream_i,
    hls::stream<ap_uint<SMER_SIZE>>& stream_o,
    ap_uint<SMER_SIZE>& debug_first_smer,
    ap_uint<SMER_SIZE>& debug_first_hash
) {
#pragma HLS INLINE off
    ap_uint<SMER_SIZE> fwd = 0;
    ap_uint<SMER_SIZE> rev = 0;

    bool first_smer_done = false;

    for (int i = 0; i < SMER-1; i++) {
#pragma HLS PIPELINE II=1
        uint8_t c = stream_i.read();
        if (c == 0) { stream_o.write(0); return; }
        ap_uint<2> b = nucl_encode(c);
        fwd <<= 2; fwd(1,0) = b;
        rev >>= 2; rev(SMER_SIZE-1, SMER_SIZE-2) = (0x2 ^ b);
    }

    while (true) {
#pragma HLS PIPELINE II=1
        uint8_t c = stream_i.read();
        if (c == 0) break;

        ap_uint<2> b = nucl_encode(c);
        fwd <<= 2; fwd(1,0) = b;
        rev = (rev >> 2) | ((ap_uint<SMER_SIZE>)(0x2 ^ b) << (SMER_SIZE-2));

        ap_uint<SMER_SIZE> canon = min_v1(fwd, rev);
        ap_uint<64> hash_val = hash_u64(canon);

        // Sauvegarde le premier s-mer et hash pour debug
        if (!first_smer_done) {
            debug_first_smer = canon;
            debug_first_hash = hash_val;
            first_smer_done = true;
        }

        stream_o.write(hash_val);
    }

    stream_o.write(0);
}

// ============================================================
// THREAD DEDUP
// ============================================================
void thread_dedup(
    hls::stream<ap_uint<SMER_SIZE>>& stream_i,
    hls::stream<ap_uint<SMER_SIZE>>& stream_o
) {
#pragma HLS INLINE off
    ap_uint<SMER_SIZE> buffer[WINDOW_SIZE];
#pragma HLS ARRAY_PARTITION variable=buffer complete

    for (int i = 0; i < WINDOW_SIZE-1; i++) {
#pragma HLS PIPELINE II=1
        buffer[i] = stream_i.read();
    }

    ap_uint<SMER_SIZE> last_min = (ap_uint<SMER_SIZE>)-1;
    int pos = 0;

    while (true) {
#pragma HLS PIPELINE II=1
        ap_uint<SMER_SIZE> v = stream_i.read();
        if (v == 0) break;

        buffer[(WINDOW_SIZE-1 + pos) % WINDOW_SIZE] = v;

        ap_uint<SMER_SIZE> m = buffer[pos];
        for (int i = 1; i < WINDOW_SIZE; i++) {
#pragma HLS UNROLL
            int idx = (pos + i) % WINDOW_SIZE;
            if (buffer[idx] < m) m = buffer[idx];
        }

        pos = (pos + 1) % WINDOW_SIZE;

        if (m != last_min) {
            last_min = m;
            stream_o.write(m);
        }
    }

    stream_o.write(0);
}

// ============================================================
// THREAD STORE
// ============================================================
void thread_store(
    hls::stream<ap_uint<SMER_SIZE>>& stream_i,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nElements
) {
#pragma HLS INLINE off
    int cnt = 0;
    while (true) {
#pragma HLS PIPELINE II=1
        ap_uint<SMER_SIZE> v = stream_i.read();
        if (v == 0) break;
        tab_hash[cnt++] = v;
    }
    *nElements = cnt;
}

// ============================================================
// TOP FUNCTION
// ============================================================
void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs,

    ap_uint<SMER_SIZE>& debug_first_smer,
    ap_uint<SMER_SIZE>& debug_first_hash
) {
#pragma HLS INTERFACE mode=m_axi     port=packed_sequence
#pragma HLS INTERFACE mode=m_axi     port=tab_hash
#pragma HLS INTERFACE mode=s_axilite port=n
#pragma HLS INTERFACE mode=s_axilite port=nMinizrs
#pragma HLS INTERFACE mode=s_axilite port=debug_first_smer
#pragma HLS INTERFACE mode=s_axilite port=debug_first_hash
#pragma HLS INTERFACE mode=s_axilite port=return

#pragma HLS DATAFLOW
    hls::stream<uint8_t, DATA_DEPTH> fifo_1;
    hls::stream<ap_uint<SMER_SIZE>, DATA_DEPTH> fifo_2;
    hls::stream<ap_uint<SMER_SIZE>, DATA_DEPTH> fifo_3;

#pragma HLS STREAM variable=fifo_1 depth=1024
#pragma HLS STREAM variable=fifo_2 depth=1024
#pragma HLS STREAM variable=fifo_3 depth=1024

    // Pipeline avec debug
    thread_reader(packed_sequence, n ,fifo_1);
    thread_smer(fifo_1, fifo_2, debug_first_smer, debug_first_hash);
    thread_dedup(fifo_2, fifo_3);
    thread_store(fifo_3, tab_hash, nMinizrs);
}
