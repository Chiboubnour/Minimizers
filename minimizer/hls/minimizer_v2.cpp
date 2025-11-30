#include <ap_int.h>
#include <hls_stream.h>

#define SMER_SIZE         56
#define DATA_DEPTH        1024
#define SMERS_PER_CYCLE   8
#define DEDUP_WINDOW_SIZE 16

inline ap_uint<2> nucl_encode(char nucl) {
    switch (nucl) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default : return 0;
    }
}

inline ap_uint<64> min_v1(const ap_uint<64> a, const ap_uint<64> b) {
    return (a < b) ? a : b;
}

inline ap_uint<64> mask_right(int numbits) {
    return (numbits >= 64) ? ~0ULL : ((1ULL << numbits) - 1ULL);
}

inline ap_uint<64> hash_u64(ap_uint<64> key, ap_uint<64> mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

void thread_reader_v2(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n_bases,
    hls::stream< ap_uint<24> >& stream_o
) {
    const int n_words = (int)((n_bases + 7) / 8);
    // Boucle principale : lecture mot par mot (8 bases par mot 64 bits)
    for (int i = 0; i < n_words; ++i) {
    #pragma HLS PIPELINE II=1
        const ap_uint<64> word_8b = packed_sequence[i];
        ap_uint<24> word_3b = 0;

        bool all_valid = true;

        // Encodage des 8 bases en 3 bits (2 bits code + 1 bit validité)
        for (int j = 0; j < 8; ++j) {
        #pragma HLS UNROLL
            const int idx = i * 8 + j;
            ap_uint<1> valid = (idx < (int)n_bases) ? 1 : 0;

            ap_uint<8> c = 0;
            if (valid) {
                c = word_8b.range(8*(j+1)-1, 8*j);
            }

            ap_uint<2> enc = 0;
            if (valid) {
                enc = nucl_encode((char)c);
            }

            word_3b.range(3*j+1, 3*j) = enc;
            word_3b[3*j+2] = valid;

            all_valid &= (bool)valid;
        }

        stream_o.write(word_3b);

        if (!all_valid)
            break;
    }

    if ((n_bases & 7) == 0)
        stream_o.write((ap_uint<24>)0);
}



void thread_reader_pack(
    hls::stream<ap_uint<24>>& stream_i,
    hls::stream<ap_uint<SMERS_PER_CYCLE*3>>& stream_o
) {
    ap_uint<SMERS_PER_CYCLE*3> pack = 0;
    int count = 0;

    while (true) {
        #pragma HLS PIPELINE II=1
        ap_uint<24> val = stream_i.read();

        if (val == 0) {
            if (count != 0)
                stream_o.write(pack);  
            stream_o.write(0);
            break;
        }

        for (int i = 0; i < 8; i++) {
            pack.range((count+1)*3 - 1, count*3) = val.range(3*i+1, 3*i);
            count++;

            if (count == SMERS_PER_CYCLE) {
                stream_o.write(pack);
                count = 0;
                pack = 0;
            }
        }
    }
}

void thread_smer_v2(
    hls::stream<ap_uint<SMERS_PER_CYCLE*3>>& stream_i,
    ap_uint<64> n_bases,
    hls::stream<ap_uint<SMER_SIZE>>& stream_o
) {
    constexpr int smer_bases = SMER_SIZE / 2;
    const ap_uint<64> HASH_MASK = mask_right(SMER_SIZE);

    ap_uint<SMER_SIZE> current_smer = 0;
    ap_uint<SMER_SIZE> cur_inv_smer = 0;
    ap_uint<64> base_count = 0; // number of bases processed so far

    while (true) {
    #pragma HLS PIPELINE II=1
        ap_uint<SMERS_PER_CYCLE*3> word = stream_i.read();
        if (word == 0) {
            stream_o.write(0);
            break;
        }

        for (int i = 0; i < SMERS_PER_CYCLE; ++i) {
        #pragma HLS UNROLL
            const ap_uint<2> c_nucl = word.range(3*i+1, 3*i);

            current_smer = (current_smer << 2);
            current_smer(1,0) = c_nucl;

            cur_inv_smer = (cur_inv_smer >> 2) |
                ((ap_uint<SMER_SIZE>)(((ap_uint<2>)0x2) ^ c_nucl) << (SMER_SIZE - 2));

            base_count++;

            if (base_count >= ( ap_uint<64>)smer_bases) {
                ap_uint<SMER_SIZE> vmin = min_v1(current_smer, cur_inv_smer);
                ap_uint<SMER_SIZE> vhash = hash_u64((ap_uint<64>)vmin, (ap_uint<64>)HASH_MASK);
                stream_o.write(vhash);
            }
        }
    }
}

void thread_smer_pack(
    hls::stream<ap_uint<SMER_SIZE>>& stream_i,
    hls::stream<ap_uint<SMER_SIZE*SMERS_PER_CYCLE>>& stream_o
) {
    ap_uint<SMER_SIZE*SMERS_PER_CYCLE> pack = 0;
    int count = 0;

    while (true) {
        #pragma HLS PIPELINE II=1

        ap_uint<SMER_SIZE> val = stream_i.read();

        // Fin de flux
        if (val == 0) {
            if (count != 0)
                stream_o.write(pack);
            stream_o.write(0);
            break;
        }

        pack.range((count+1)*SMER_SIZE - 1, count*SMER_SIZE) = val;
        count++;

        if (count == SMERS_PER_CYCLE) {
            stream_o.write(pack);
            count = 0;
            pack = 0;
        }
    }
}

void thread_dedup_v2(
    hls::stream<ap_uint<SMER_SIZE*SMERS_PER_CYCLE>>& stream_i,
    hls::stream<ap_uint<SMER_SIZE>>& stream_o
) {
    ap_uint<SMER_SIZE> buffer[DEDUP_WINDOW_SIZE];
    #pragma HLS ARRAY_PARTITION variable=buffer complete dim=1

    for (int i = 0; i < DEDUP_WINDOW_SIZE; i++)
        buffer[i] = (ap_uint<SMER_SIZE>)(-1);

    ap_uint<SMER_SIZE> lastOutput = (ap_uint<SMER_SIZE>)(-1);

    while (true) {
    #pragma HLS PIPELINE II=1

        ap_uint<SMER_SIZE*SMERS_PER_CYCLE> pack = stream_i.read();
        if (pack == 0) {
            stream_o.write(0);
            break;
        }

        for (int i = 0; i < SMERS_PER_CYCLE; i++) {
        #pragma HLS UNROLL

            ap_uint<SMER_SIZE> vhash = pack.range((i+1)*SMER_SIZE - 1, i*SMER_SIZE);

            ap_uint<SMER_SIZE> minz = vhash;
            for (int p = 0; p < DEDUP_WINDOW_SIZE; p++) {
                #pragma HLS UNROLL
                if (buffer[p] < minz)
                    minz = buffer[p];
            }
            for (int p = 0; p < DEDUP_WINDOW_SIZE-1; p++) {
                #pragma HLS UNROLL
                buffer[p] = buffer[p+1];
            }
            buffer[DEDUP_WINDOW_SIZE-1] = vhash;

            if (minz != lastOutput) {
                stream_o.write(minz);
                lastOutput = minz;
            }
        }
    }
}


void thread_store_v2(
    hls::stream<ap_uint<SMER_SIZE>>& stream_i,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nElements
){
    int cnt=0;
    while(true){
    #pragma HLS PIPELINE II=1
        ap_uint<SMER_SIZE> val=stream_i.read();
        if(val==0) break;
        tab_hash[cnt++]=val;
    }
    *nElements=cnt;
}

void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs
){
    #pragma HLS INTERFACE m_axi port=packed_sequence      offset=slave bundle=gmem0
    #pragma HLS INTERFACE m_axi port=tab_hash             offset=slave bundle=gmem1
    #pragma HLS INTERFACE s_axilite port=n                bundle=control
    #pragma HLS INTERFACE s_axilite port=nMinizrs         bundle=control
    #pragma HLS INTERFACE s_axilite port=return           bundle=control
	#pragma HLS DATAFLOW

    hls::stream< ap_uint<24>, DATA_DEPTH > fifo_1;
    hls::stream< ap_uint<SMERS_PER_CYCLE*3>, DATA_DEPTH > fifo_2;
    hls::stream< ap_uint<SMER_SIZE>, DATA_DEPTH > fifo_3;
    hls::stream< ap_uint<SMER_SIZE*SMERS_PER_CYCLE>, DATA_DEPTH > fifo_4;
    hls::stream< ap_uint<SMER_SIZE>, DATA_DEPTH > fifo_5;

    hls::stream<ap_uint<64>, DATA_DEPTH> fifo_in;

    thread_reader_v2(packed_sequence, n, fifo_1);

    thread_reader_pack(fifo_1, fifo_2);

    thread_smer_v2(fifo_2, n, fifo_3);

    thread_smer_pack(fifo_3, fifo_4);

    thread_dedup_v2(fifo_4, fifo_5);

    thread_store_v2(fifo_5, tab_hash, nMinizrs);
}
