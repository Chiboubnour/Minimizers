#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

constexpr int SMER      = 19;
constexpr int WINDOW    = 16;

inline ap_uint<2> nucl_encode(uint8_t c)
{
#pragma HLS INLINE
    return (c >> 1) & 0x3;
}

template<int width>
inline ap_uint<64>min_v1(const ap_uint<width> a, const ap_uint<width> b)
{
    #pragma HLS INLINE
    return (a < b) ? a : b;
}

template <int resu_size>
inline ap_uint<resu_size> hash_u64(ap_uint<64> key) {
    key = (~key + (key << 21)) /*& mask*/; // key = (key << 21) - key - 1;
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8)) /*& mask*/; // key * 265
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4)) /*& mask*/; // key * 21
    key ^= key >> 28;
    key = (key + (key << 31)) /*& mask*/;
    return key;
}

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
        ap_uint<8> valid;
        if (w == n_words - 1 && (n_bases & 7)) {
            valid = n_bases & 7;
        } else {
            valid = 8;
        }

        base_stream.write(word);
        base_valid.write(valid);
    }
    base_stream.write(0);
    base_valid.write(0);
}

static void adapter_hls(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<64>>& out_stream,
    hls::stream<ap_uint<8>>&  out_valid
){
#pragma HLS INLINE off

    ap_uint<128> buffer = 0;
    ap_uint<8> valid_count = 0;


    while (true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> in_word  = base_stream.read();
        ap_uint<8>  in_valid = base_valid.read();
        if (in_valid == 0){

            if(valid_count > 0){
                ap_uint<64> out_word = buffer.range(63,0);
                out_stream.write(out_word);
                out_valid.write(valid_count);
            }


            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        buffer.range(valid_count*8 + 63, valid_count*8) = in_word;
        valid_count += in_valid;

        if(valid_count >= 8){

            ap_uint<64> out_word = buffer.range(63,0);
            out_stream.write(out_word);
            out_valid.write(8);

            buffer >>= 64;
            valid_count -= 8;
        }
    }
}

template <int S_BASE>
void smer_thread_hls(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<16*S_BASE>>& out_stream,
    hls::stream<ap_uint<8>>&  out_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * S_BASE;

    ap_uint<S_BITS> current_smer = 0;
    ap_uint<S_BITS> cur_inv_smer = 0;
    int base_count = 0;

    ap_uint<8*S_BITS> buffer = 0;
    ap_uint<8> count_smers = 0;

    while (true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if(valid == 0){
            if(count_smers > 0){
                out_stream.write(buffer);
                out_valid.write(count_smers);
            }

            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        // ---- boucle FIXE 8 ----
        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            if(i < valid){

                uint8_t nucl = (word >> (8*i)) & 0xFF;
                ap_uint<2> c_nucl = (nucl >> 1) & 0x3;

                // rolling smer
                current_smer <<= 2;
                current_smer(1,0) = c_nucl;

                cur_inv_smer >>= 2;
                cur_inv_smer(S_BITS-1, S_BITS-2) = 0x2 ^ c_nucl;

                base_count++;

                if(base_count >= S_BASE){

                    ap_uint<S_BITS> vmin  =
                        (current_smer < cur_inv_smer)
                        ? current_smer
                        : cur_inv_smer;

                    ap_uint<S_BITS> vhash = hash_u64<S_BITS>((ap_uint<64>)vmin);

                    buffer.range(
                        (count_smers+1)*S_BITS-1,
                        count_smers*S_BITS
                    ) = vhash;

                    count_smers++;

                    if(count_smers == 8){
                        out_stream.write(buffer);
                        out_valid.write(8);
                        buffer = 0;
                        count_smers = 0;
                    }
                }
            }
        }
    }

    if(count_smers > 0){
        out_stream.write(buffer);
        out_valid.write(count_smers);
    }
}


template<int S_BITS>
static void adapter_smer_hls(
    hls::stream<ap_uint<8*S_BITS>>& base_stream,
    hls::stream<ap_uint<8>>&        base_valid,
    hls::stream<ap_uint<8*S_BITS>>& out_stream,
    hls::stream<ap_uint<8>>&        out_valid
)
{
#pragma HLS INLINE off

    ap_uint<16*S_BITS> buffer = 0;
    ap_uint<4> valid_count = 0;

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto in_word  = base_stream.read();
        auto in_valid = base_valid.read();

        // ===== EOS =====
        if(in_valid == 0){

            if(valid_count > 0){
                out_stream.write(buffer.range(8*S_BITS-1,0));
                out_valid.write(valid_count);
            }

            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        buffer |= ((ap_uint<16*S_BITS>)in_word) << (valid_count*S_BITS);
        valid_count += in_valid;

        if(valid_count >= 8){

            out_stream.write(buffer.range(8*S_BITS-1,0));
            out_valid.write(8);

            buffer >>= 8*S_BITS;
            valid_count -= 8;
        }
    }
}

template<int WINDOW,int SMER>
void thread_dedup_v8(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> dq_val[WINDOW];
    ap_uint<4>    dq_pos[WINDOW];
#pragma HLS ARRAY_PARTITION variable=dq_val complete
#pragma HLS ARRAY_PARTITION variable=dq_pos complete

    ap_uint<4> head = 0;
    ap_uint<4> tail = 0;

    ap_uint<4> cur_pos = 0;
    ap_uint<SMER> last = ~ap_uint<SMER>(0);

    for(int p=0;p < 2;p++){
#pragma HLS PIPELINE II=1

        auto packed = in.read();
        auto valid  = valid_i.read();
        if(valid==0){
            out.write(0);
            valid_o.write(0);
            return;
        }

        for(int i=0;i<8;i++){
#pragma HLS PIPELINE II=1

            ap_uint<SMER> v =
                packed.range((i+1)*SMER-1, i*SMER);

            dq_val[tail] = v;
            dq_pos[tail] = cur_pos;
            tail++;
            cur_pos++;
        }
    }

    for(int i=0;i<WINDOW;i++){
        ap_uint<SMER> v = dq_val[i];

        while(head != tail && dq_val[(tail-1)&15] > v){
            tail--;
        }

        dq_val[tail] = v;
        dq_pos[tail] = i;
        tail++;
    }

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto packed = in.read();
        auto valid  = valid_i.read();

        if(valid==0){
            out.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<8*SMER> pout = 0;
        ap_uint<4> ocnt = 0;

        for(int i=0;i<8;i++){

            if(i>=valid) break;

            ap_uint<SMER> v =
                packed.range((i+1)*SMER-1, i*SMER);

            if(dq_pos[head] == (cur_pos - WINDOW)){
                head++;
            }

            while(head != tail &&
                  dq_val[(tail-1)&15] > v){
                tail--;
            }

            dq_val[tail] = v;
            dq_pos[tail] = cur_pos;
            tail++;

            ap_uint<SMER> m = dq_val[head];

            if(m != last){
                last = m;
                pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                ocnt++;
            }

            cur_pos++;
        }

        if(ocnt>0){
            out.write(pout);
            valid_o.write(ocnt);
        }
    }
}


template<int SMER>
void thread_store_burst_v8(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>&       valid_i,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nElem
){
#pragma HLS INLINE off

    ap_uint<64> cnt = 0;

    ap_uint<64> burst_buf[64];
#pragma HLS ARRAY_PARTITION variable=burst_buf complete

    ap_uint<6> bcnt = 0;   // 0..63

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto packed = in.read();
        auto valid  = valid_i.read();

        if (valid == 0) break;

        for(int i=0;i<8;i++){

            if(i < valid){

                ap_uint<SMER> v =
                    packed.range((i+1)*SMER-1, i*SMER);

                burst_buf[bcnt] = v;
                bcnt++;

                if(bcnt == 64){

                    for(int j=0;j<64;j++){
#pragma HLS PIPELINE II=1
                        tab_hash[cnt++] = burst_buf[j];
                    }

                    bcnt = 0;
                }
            }
        }
    }

    for(int j=0;j<bcnt;j++){
#pragma HLS PIPELINE II=1
        tab_hash[cnt++] = burst_buf[j];
    }

    *nElem = cnt;
}


void minimizer(
    const ap_uint<64>* packed_sequence,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nMinizrs,
    ap_uint<64>  n_bases
)
{
#pragma HLS INTERFACE m_axi     port=packed_sequence offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi     port=tab_hash        offset=slave bundle=gmem1

#pragma HLS INTERFACE s_axilite port=nMinizrs
#pragma HLS INTERFACE s_axilite port=n_bases
#pragma HLS INTERFACE s_axilite port=return

#pragma HLS DATAFLOW

    constexpr int S_BITS = 2 * SMER;

    hls::stream<ap_uint<64>> st_base_raw;
    hls::stream<ap_uint<8>>  st_base_valid;

    hls::stream<ap_uint<64>> st_base_aligned;
    hls::stream<ap_uint<8>>  st_base_aligned_valid;

    hls::stream<ap_uint<8*S_BITS>> st_smer_packed;
    hls::stream<ap_uint<8>>        st_smer_packed_valid;

    hls::stream<ap_uint<8*S_BITS>> st_smer_aligned;
    hls::stream<ap_uint<8>>        st_smer_aligned_valid;

    hls::stream<ap_uint<8*S_BITS>> st_minimizers;
    hls::stream<ap_uint<8>>        st_minimizers_valid;

#pragma HLS STREAM variable=st_base_raw           depth=256
#pragma HLS STREAM variable=st_base_valid         depth=256
#pragma HLS STREAM variable=st_base_aligned       depth=256
#pragma HLS STREAM variable=st_base_aligned_valid depth=256
#pragma HLS STREAM variable=st_smer_packed        depth=256
#pragma HLS STREAM variable=st_smer_packed_valid  depth=256
#pragma HLS STREAM variable=st_smer_aligned       depth=256
#pragma HLS STREAM variable=st_smer_aligned_valid depth=256
#pragma HLS STREAM variable=st_minimizers         depth=256
#pragma HLS STREAM variable=st_minimizers_valid   depth=256


    reader_hls(
        packed_sequence,
        n_bases,
        st_base_raw,
        st_base_valid
    );

    adapter_hls(
        st_base_raw,
        st_base_valid,
        st_base_aligned,
        st_base_aligned_valid
    );

    smer_thread_hls<SMER>(
        st_base_aligned,
        st_base_aligned_valid,
        st_smer_packed,
        st_smer_packed_valid
    );

    adapter_smer_hls<S_BITS>(
        st_smer_packed,
        st_smer_packed_valid,
        st_smer_aligned,
        st_smer_aligned_valid
    );

    thread_dedup_v8<WINDOW, S_BITS>(
        st_smer_aligned,
        st_smer_aligned_valid,
        st_minimizers,
        st_minimizers_valid
    );

    thread_store_burst_v8<S_BITS>(
        st_minimizers,
        st_minimizers_valid,
        tab_hash,
        nMinizrs
    );
}
