#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

constexpr int SMER      = 19;
constexpr int WINDOW    = 16;

// ================= helpers =================
inline ap_uint<2> nucl_encode(uint8_t c)
{
#pragma HLS INLINE
    return (c >> 1) & 0x3;
}

inline ap_uint<64> hash_u64(ap_uint<2*SMER> key)
{
#pragma HLS INLINE
    ap_uint<64> k64 = key;

    k64 = ~k64 + (k64 << 21);
    k64 ^= k64 >> 24;
    k64 = (k64 + (k64 << 3)) + (k64 << 8);
    k64 ^= k64 >> 14;
    k64 = (k64 + (k64 << 2)) + (k64 << 4);
    k64 ^= k64 >> 28;
    k64 = k64 + (k64 << 31);

    return k64;
}

// ================= reader =================
static void reader_hls(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n_bases,
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid
) {
#pragma HLS INLINE off
    ap_uint<64> n_words = (n_bases + 7) >> 3;

    for (ap_uint<64> w = 0; w < n_words; w++) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> word = packed_sequence[w];
        ap_uint<64> bases = word;
        ap_uint<8>  valid = 0xFF;

        if (w == n_words-1 && (n_bases & 0x7)) {
            int rem = n_bases & 0x7;
            valid = (1 << rem) - 1;
        }

        base_stream.write(bases);
        base_valid.write(valid);
    }

    base_stream.write(0);
    base_valid.write(0);
}

// ================= thread_smer_generator_v2 =================
template<int SMER = 19>
void thread_smer_generator_v2(
    hls::stream<ap_uint<64>>& base_i,
    hls::stream<ap_uint<8>>&  valid_i,
    hls::stream<ap_uint<8*SMER>>& smer_o,
    hls::stream<ap_uint<8>>&  valid_o
) {
#pragma HLS INLINE off

    constexpr int INNER_BITS = 2 * SMER;

    ap_uint<INNER_BITS> current_smer = 0;
    ap_uint<INNER_BITS> cur_inv_smer = 0;
    ap_uint<6> init_count = 0;
    bool initialized = false;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        ap_uint<64> bases = base_i.read();
        ap_uint<8>  valid = valid_i.read();

        if (valid == 0) {
            smer_o.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<8*SMER> packed_smers = 0;
        ap_uint<8> smer_count = 0;

        for (int i = 0; i < 8; i++) {
#pragma HLS UNROLL
            if (!valid[i]) continue;

            ap_uint<8> nucl = bases.range(8*i + 7, 8*i);
            if (nucl == 0) continue;

            ap_uint<2> c_nucl = (nucl >> 1) & 0x03;

            current_smer <<= 2;
            cur_inv_smer >>= 2;
            current_smer(1,0) = c_nucl;
            cur_inv_smer(INNER_BITS-1, INNER_BITS-2) = 0x2 ^ c_nucl;

            if (!initialized) {
                init_count++;
                if (init_count == SMER-1)
                    initialized = true;
                continue;
            }

            ap_uint<INNER_BITS> vmin = (current_smer < cur_inv_smer) ? current_smer : cur_inv_smer;
            ap_uint<SMER> vhash = hash_u64(vmin).range(SMER-1,0);

            packed_smers.range((smer_count+1)*SMER-1, smer_count*SMER) = vhash;
            smer_count++;
        }

        smer_o.write(packed_smers);
        valid_o.write(smer_count);
    }
}

// ================= thread_dedup_v8 =================
template <int window, int SMER>
void thread_dedup_v8(
    hls::stream<ap_uint<8*SMER>>& stream_i,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& stream_o,
    hls::stream<ap_uint<8>>& valid_o
) {
#pragma HLS INLINE off
    ap_uint<SMER> buffer[window];
#pragma HLS ARRAY_PARTITION variable=buffer complete
    int buffer_pos = 0;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        ap_uint<8*SMER> packed_smers = stream_i.read();
        ap_uint<8> valid = valid_i.read();
        if (valid == 0) {
            stream_o.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<SMER> min_vals[8];
        for (int i=0; i<8; i++) {
#pragma HLS UNROLL
            ap_uint<SMER> v = packed_smers.range((i+1)*SMER-1, i*SMER);
            buffer[buffer_pos] = v;

            ap_uint<SMER> min_val = buffer[0];
            for (int j=1; j<window; j++)
                if (buffer[j] < min_val) min_val = buffer[j];

            buffer_pos = (buffer_pos + 1) % window;
            min_vals[i] = min_val;
        }

        ap_uint<8*SMER> packed_out = 0;
        for (int i=0; i<8; i++)
            packed_out.range((i+1)*SMER-1, i*SMER) = min_vals[i];

        stream_o.write(packed_out);
        valid_o.write(8);
    }
}

// ================= thread_store_burst_v8 =================
template<int STREAM_BITS>
void thread_store_burst_v8(
    hls::stream<ap_uint<STREAM_BITS>>& dedup_i,
    hls::stream<ap_uint<8>>& valid_i,
    ap_uint<512>* tab_hash,
    ap_uint<64>* nElements
)
{
#pragma HLS INLINE off
    ap_uint<64> cnt = 0;
    constexpr int LOCAL_SMER = STREAM_BITS / 8;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<STREAM_BITS> packed = dedup_i.read();
        ap_uint<8> valid = valid_i.read();
        if (valid == 0) break;

        ap_uint<512> word = 0;
        for (int i = 0; i < 8; i++) {
#pragma HLS UNROLL
            ap_uint<LOCAL_SMER> mini = packed.range((i+1)*LOCAL_SMER-1, i*LOCAL_SMER);
            word.range((i+1)*LOCAL_SMER-1, i*LOCAL_SMER) = mini;
        }
        tab_hash[cnt >> 3] = word;
        cnt += 8;
    }

    *nElements = cnt;
}


void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<512>* tab_hash,       // stocke 8 s-mers de 64 bits par mot
    ap_uint<64>* nMinizrs,
    ap_uint<64>  n_bases
)
{
#pragma HLS INTERFACE m_axi     port=packed_sequence   depth=256
#pragma HLS INTERFACE m_axi     port=tab_hash          depth=256
#pragma HLS INTERFACE s_axilite port=nMinizrs
#pragma HLS INTERFACE s_axilite port=n_bases
#pragma HLS INTERFACE s_axilite port=return

    // Streams internes
    hls::stream<ap_uint<64>> base_stream("base_stream");
    hls::stream<ap_uint<8>>  base_valid("base_valid");
    hls::stream<ap_uint<8*SMER>> smer_stream("smer_stream");
    hls::stream<ap_uint<8>> smer_valid("smer_valid");
    hls::stream<ap_uint<8*SMER>> dedup_stream("dedup_stream");
    hls::stream<ap_uint<8>> dedup_valid("dedup_valid");

#pragma HLS STREAM variable=base_stream  depth=256
#pragma HLS STREAM variable=base_valid  depth=256
#pragma HLS STREAM variable=smer_stream  depth=256
#pragma HLS STREAM variable=smer_valid  depth=256
#pragma HLS STREAM variable=dedup_stream depth=256
#pragma HLS STREAM variable=dedup_valid depth=256

#pragma HLS DATAFLOW
    reader_hls(packed_sequence, n_bases, base_stream, base_valid);
    thread_smer_generator_v2<SMER>(base_stream, base_valid, smer_stream, smer_valid);
    thread_dedup_v8<WINDOW, SMER>(smer_stream, smer_valid, dedup_stream, dedup_valid);
    thread_store_burst_v8<8*SMER>(dedup_stream, dedup_valid, tab_hash, nMinizrs);
}
