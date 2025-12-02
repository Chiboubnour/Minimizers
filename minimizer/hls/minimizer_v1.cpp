#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

#define SMER    19
#define WINDOW  16
#define SMER_SIZE (2*SMER)

template<int W>
static inline ap_uint<W> min_v1(const ap_uint<W> a, const ap_uint<W> b) {
#pragma HLS INLINE
    return (a < b) ? a : b;
}

template<int W>
static inline ap_uint<W> hash_u64(ap_uint<W> key) {
#pragma HLS INLINE
    key = (~key + (key << 21));
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >> 28);
    key = (key + (key << 31));
    return key;
}

static inline ap_uint<2> nucl_encode(const ap_uint<8> nucl) {
#pragma HLS INLINE
    return (nucl >> 1) & 0x3;
}

template<int SMER_SZ = SMER>
void thread_reader(const ap_uint<64>* packed_sequence, hls::stream<uint8_t>& out) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
    bool stop = false;
    // on lit jusqu'à 64 mots (comme la référence)
    for (int i = 0; i < 64 && !stop; ++i) {
        ap_uint<64> word = packed_sequence[i];
        for (int j = 0; j < 8; ++j) {
            uint8_t c = (uint8_t)(word & 0xFF);
            if (c == 0) {
                stop = true;
                break;
            }
            out.write(c);
            word >>= 8;
        }
    }
    // End marker
    out.write((uint8_t)0);
}

template <int SMER_BITS>
void thread_smer(hls::stream<uint8_t>& in, hls::stream<ap_uint<SMER_BITS>>& out) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1

    ap_uint<SMER_BITS> current = 0;
    ap_uint<SMER_BITS> cur_inv = 0;
    const int smer = SMER_BITS / 2;

    // Build first s-mer (smer-1 bases)
    for (int i = 0; i < smer - 1; ++i) {
        uint8_t nu = in.read();
        if (nu == 0) {
            // insufficient length -> terminate with sentinel
            out.write((ap_uint<SMER_BITS>)0);
            return;
        }
        ap_uint<2> c = nucl_encode(nu);
        current <<= 2;
        current.range(1,0) = c;
        cur_inv >>= 2;
        cur_inv.range(SMER_BITS-1, SMER_BITS-2) = (0x2 ^ c);
    }

    // Main loop: read until 0 sentinel
    while (true) {
        uint8_t nu = in.read();
        if (nu == 0) break;
        ap_uint<2> c = nucl_encode(nu);

        current <<= 2;
        current.range(1,0) = c;

        cur_inv = (cur_inv >> 2) | ((ap_uint<SMER_BITS>)((0x2 ^ c)) << (SMER_BITS-2));

        const ap_uint<SMER_BITS> vmin = min_v1<SMER_BITS>(current, cur_inv);
        const ap_uint<SMER_BITS> vhash = hash_u64<SMER_BITS>(vmin);

        out.write(vhash);
    }

    // End marker for s-mer stream
    out.write((ap_uint<SMER_BITS>)0);
}

template <int WINDOW_SZ, int SMER_BITS>
void thread_dedup(
    hls::stream<ap_uint<SMER_BITS>>& in,
    hls::stream<ap_uint<SMER_BITS>>& out
)
                  {

#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_PARTITION variable=buffer complete

    ap_uint<SMER_BITS> buffer[WINDOW_SZ];
    ap_uint<SMER_BITS> last_mini = (ap_uint<SMER_BITS>)(-1);
    int pos = 0;

    for (int i = 0; i < WINDOW_SZ - 1; ++i) {
        ap_uint<SMER_BITS> v = in.read();
        buffer[i] = v;
    }

    while (true) {
        ap_uint<SMER_BITS> v = in.read();
        if (v == 0) break;

        buffer[(WINDOW_SZ - 1 + pos) % WINDOW_SZ] = v;

        ap_uint<SMER_BITS> min_val = buffer[pos];
        for (int k = 1; k < WINDOW_SZ; ++k) {
            int idx = (pos + k) % WINDOW_SZ;
            if (buffer[idx] < min_val) min_val = buffer[idx];
        }

        pos = (pos + 1) % WINDOW_SZ;

        if (min_val != last_mini) {
            last_mini = min_val;
            out.write(min_val);
        }
    }

    out.write((ap_uint<SMER_BITS>)0);
}

template <int SMER_BITS>
void thread_store(hls::stream<ap_uint<SMER_BITS>>& in,
                  ap_uint<64>* tab_hash,
                  ap_uint<64>* nElements) {

    int cnt = 0;
    while (true) {
        ap_uint<SMER_BITS> v = in.read();
        if (v == 0) break;
        tab_hash[cnt++] = (ap_uint<64>)v;
    }
    *nElements = (ap_uint<64>)cnt;
}


template <int SMER_N, int WINDOW_N>
void thr_minimizer(const uint64_t* packed_sequence,
                   ap_uint<64>* tab_hash,
                   ap_uint<64>* nMinizrs
                ) 
                 
                   {


#pragma HLS DATAFLOW
#pragma HLS INTERFACE m_axi port=packed_sequence offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=tab_hash offset=slave bundle=gmem
#pragma HLS INTERFACE s_axilite port=packed_sequence bundle=control
#pragma HLS INTERFACE s_axilite port=tab_hash bundle=control
#pragma HLS INTERFACE s_axilite port=nMinizrs bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    hls::stream<uint8_t> fifo1("reader_to_smer");
    hls::stream<ap_uint<2*SMER_N>> fifo2("smer_to_dedup");
    hls::stream<ap_uint<2*SMER_N>> fifo3("dedup_to_store");

#pragma HLS STREAM variable=fifo1 depth=256
#pragma HLS STREAM variable=fifo2 depth=256
#pragma HLS STREAM variable=fifo3 depth=256

    thread_reader<SMER_N>((const ap_uint<64>*)packed_sequence, fifo1);
    thread_smer<2*SMER_N>(fifo1, fifo2);
    thread_dedup<WINDOW_N, 2*SMER_N>(fifo2, fifo3);
    thread_store<2*SMER_N>(fifo3, tab_hash, nMinizrs);
}

extern "C" void minimizer_hls_top(const ap_uint<64>* packed_sequence,
                                  ap_uint<64>* tab_hash,
                                  ap_uint<64>* nMinizrs)
                                   {
    thr_minimizer<SMER, WINDOW>( (const uint64_t*)packed_sequence, tab_hash, nMinizrs );
}
