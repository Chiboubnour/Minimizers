
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

constexpr int SMER            = 21;
constexpr int WINDOW          = 11;
constexpr int Parallel_factor = 8;

inline ap_uint<2> nucl_encode(uint8_t c)
{
#pragma HLS INLINE
    return (c >> 1) & 0x3;
}

template<int width>
inline ap_uint<64> min_v1(const ap_uint<width> a, const ap_uint<width> b)
{
#pragma HLS INLINE
    return (a < b) ? a : b;
}

template <int resu_size>
inline ap_uint<resu_size> hash_u64(ap_uint<64> key) {
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return key;
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

//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
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

//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
template<int S_BASE>
void smer_generate_thread(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<Parallel_factor*2*S_BASE>>& smer_stream,
    hls::stream<ap_uint<8>>&  smer_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;
    const ap_uint<S_BITS> MASK = ap_uint<S_BITS>(-1);

    ap_uint<S_BITS> fwd_state = 0;
    ap_uint<S_BITS> rev_state = 0;
    ap_uint<32> base_cnt = 0;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if(valid == 0){
            smer_stream.write(0);
            smer_valid.write(0);
            break;
        }

        ap_uint<Parallel_factor*S_BITS> packed_out = 0;
        ap_uint<S_BITS> fwd_local = fwd_state;
        ap_uint<S_BITS> rev_local = rev_state;
        int out_idx = 0;

    BASE_LOOP:
        for(int i = 0; i < Parallel_factor; i++){

            ap_uint<8> nucl = word.range(8*i+7, 8*i);
            ap_uint<2> c    = (nucl >> 1) & 3;
            ap_uint<2> comp = 0x2 ^ c;

            fwd_local <<= 2;
            fwd_local(1,0) = c;

            rev_local >>= 2;
            rev_local(S_BITS-1, S_BITS-2) = comp;

            base_cnt++;

            if(base_cnt >= S_BASE){
                ap_uint<S_BITS> canon = (fwd_local < rev_local) ? fwd_local : rev_local;
                packed_out.range((out_idx+1)*S_BITS-1, out_idx*S_BITS) = canon;
                out_idx++;
            }
        }

        fwd_state = fwd_local & MASK;
        rev_state = rev_local;

        if(out_idx > 0){
            smer_stream.write(packed_out);
            smer_valid.write(out_idx);
        }
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
template<int S_BASE>
void hash_thread(
    hls::stream<ap_uint<Parallel_factor*2*S_BASE>>& smer_stream,
    hls::stream<ap_uint<8>>&  smer_valid,
    hls::stream<ap_uint<Parallel_factor*2*S_BASE>>& hash_stream,
    hls::stream<ap_uint<8>>&  hash_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<Parallel_factor*S_BITS> packed_smer = smer_stream.read();
        ap_uint<8> valid = smer_valid.read();

        if(valid == 0){
            hash_stream.write(0);
            hash_valid.write(0);
            break;
        }

        ap_uint<Parallel_factor*S_BITS> packed_hash = 0;

    HASH_LOOP:
        for(int i = 0; i < Parallel_factor; i++){

            if(i < valid){
                ap_uint<S_BITS> smer = packed_smer.range((i+1)*S_BITS-1, i*S_BITS);
                ap_uint<S_BITS> h    = hash_u64<S_BITS>((ap_uint<64>)smer);
                packed_hash.range((i+1)*S_BITS-1, i*S_BITS) = h;
            }
        }

        hash_stream.write(packed_hash);
        hash_valid.write(valid);
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
template<int S_BITS>
static void adapter_smer_hls(
    hls::stream<ap_uint<Parallel_factor*S_BITS>>& base_stream,
    hls::stream<ap_uint<8>>&                      base_valid,
    hls::stream<ap_uint<Parallel_factor*S_BITS>>& out_stream,
    hls::stream<ap_uint<8>>&                      out_valid
)
{
#pragma HLS INLINE off

    ap_uint<2*Parallel_factor*S_BITS> buffer = 0;
    ap_uint<4> valid_count = 0;

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto in_word  = base_stream.read();
        auto in_valid = base_valid.read();

        if(in_valid == 0){

            if(valid_count > 0){
                out_stream.write(buffer.range(Parallel_factor*S_BITS-1, 0));
                out_valid.write(valid_count);
            }

            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        buffer |= ((ap_uint<2*Parallel_factor*S_BITS>)in_word) << (valid_count*S_BITS);
        valid_count += in_valid;

        if(valid_count >= Parallel_factor){

            out_stream.write(buffer.range(Parallel_factor*S_BITS-1, 0));
            out_valid.write(Parallel_factor);

            buffer >>= Parallel_factor*S_BITS;
            valid_count -= Parallel_factor;
        }
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
template<int WINDOW, int SMER>
void thread_min_v8(
    hls::stream<ap_uint<Parallel_factor*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<Parallel_factor*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> buf[WINDOW];
#pragma HLS ARRAY_PARTITION variable=buf complete

    ap_uint<4>  pos   = 0;
    ap_uint<16> count = 0;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto packed = in.read();
        auto valid  = valid_i.read();

        if(valid == 0){
            out.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<Parallel_factor*SMER> pout = 0;
        ap_uint<4> ocnt = 0;

    LOOP_I:
        for(int i = 0; i < Parallel_factor; i++){

            bool en = (i < valid);

            ap_uint<SMER> v = packed.range((i+1)*SMER-1, i*SMER);

            if(en){
                buf[pos] = v;

                pos++;
                if(pos == WINDOW)
                    pos = 0;

                count++;
            }

            bool ready = (count >= WINDOW);

            if(en && ready){

                ap_uint<SMER> tmp[WINDOW];

                for(int j = 0; j < WINDOW; j++){
                    tmp[j] = buf[j];
                }

                for(int s = WINDOW/2; s > 0; s >>= 1){
                    for(int j = 0; j < s; j++){
                        tmp[j] = (tmp[2*j] < tmp[2*j+1]) ?
                                  tmp[2*j] : tmp[2*j+1];
                    }
                }

                ap_uint<SMER> m = tmp[0];

                pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                ocnt++;
            }
        }

        if(ocnt != 0){
            out.write(pout);
            valid_o.write(ocnt);
        }
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

template<int WINDOW, int SMER>
void thread_dedup_v8(
    hls::stream<ap_uint<Parallel_factor*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<Parallel_factor*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> last = ~ap_uint<SMER>(0);

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto packed = in.read();
        auto valid  = valid_i.read();

        if(valid == 0){
            out.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<Parallel_factor*SMER> pout = 0;
        ap_uint<4> ocnt = 0;

        for(int i = 0; i < Parallel_factor; i++){

            ap_uint<SMER> v = packed.range((i+1)*SMER-1, i*SMER);

            if( (v != last) && (i < valid) ){
                last = v;
                pout.range((ocnt+1)*SMER-1, ocnt*SMER) = v;
                ocnt++;
            }
        }

        if(ocnt > 0){
            out.write(pout);
            valid_o.write(ocnt);
        }
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
template<int SMER>
void thread_store_burst_v8(
    hls::stream<ap_uint<Parallel_factor*SMER>>& in,
    hls::stream<ap_uint<8>>&                    valid_i,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nElem
){
#pragma HLS INLINE off

    ap_uint<64> cnt = 0;

    ap_uint<64> burst_buf[64];
#pragma HLS ARRAY_PARTITION variable=burst_buf complete

    ap_uint<6> bcnt = 0;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        auto packed = in.read();
        auto valid  = valid_i.read();

        if (valid == 0) break;

        for(int i = 0; i < Parallel_factor; i++){

            if(i < valid){

                ap_uint<SMER> v = packed.range((i+1)*SMER-1, i*SMER);

                burst_buf[bcnt] = v;
                bcnt++;

                if(bcnt == 64){

                    for(int j = 0; j < 64; j++){
#pragma HLS PIPELINE II=1
                        tab_hash[cnt++] = burst_buf[j];
                    }

                    bcnt = 0;
                }
            }
        }
    }

    for(int j = 0; j < bcnt; j++){
#pragma HLS PIPELINE II=1
        tab_hash[cnt++] = burst_buf[j];
    }

    *nElem = cnt;
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

    hls::stream<ap_uint<Parallel_factor*S_BITS>> st_smer_raw;
    hls::stream<ap_uint<8>>                      st_smer_raw_valid;

    hls::stream<ap_uint<Parallel_factor*S_BITS>> st_smer_packed;
    hls::stream<ap_uint<8>>                      st_smer_packed_valid;

    hls::stream<ap_uint<Parallel_factor*S_BITS>> st_smer_aligned;
    hls::stream<ap_uint<8>>                      st_smer_aligned_valid;

    hls::stream<ap_uint<Parallel_factor*S_BITS>> st_mins;
    hls::stream<ap_uint<8>>                      st_mins_valid;

    hls::stream<ap_uint<Parallel_factor*S_BITS>> st_minimizers;
    hls::stream<ap_uint<8>>                      st_minimizers_valid;

#pragma HLS STREAM variable=st_base_raw           depth=2
#pragma HLS STREAM variable=st_base_valid         depth=2
#pragma HLS STREAM variable=st_base_aligned       depth=2
#pragma HLS STREAM variable=st_base_aligned_valid depth=2

#pragma HLS STREAM variable=st_smer_raw           depth=2
#pragma HLS STREAM variable=st_smer_raw_valid     depth=2

#pragma HLS STREAM variable=st_smer_packed        depth=2
#pragma HLS STREAM variable=st_smer_packed_valid  depth=2

#pragma HLS STREAM variable=st_smer_aligned       depth=2
#pragma HLS STREAM variable=st_smer_aligned_valid depth=2

#pragma HLS STREAM variable=st_mins               depth=2
#pragma HLS STREAM variable=st_mins_valid         depth=2

#pragma HLS STREAM variable=st_minimizers         depth=2
#pragma HLS STREAM variable=st_minimizers_valid   depth=2

    reader_hls(packed_sequence, n_bases, st_base_raw, st_base_valid);

    adapter_hls(st_base_raw, st_base_valid, st_base_aligned, st_base_aligned_valid);

    smer_generate_thread<SMER>(st_base_aligned, st_base_aligned_valid, st_smer_raw, st_smer_raw_valid);

    hash_thread<SMER>(st_smer_raw, st_smer_raw_valid, st_smer_packed, st_smer_packed_valid);

    adapter_smer_hls<S_BITS>(st_smer_packed, st_smer_packed_valid, st_smer_aligned, st_smer_aligned_valid);

    thread_min_v8<WINDOW, S_BITS>(st_smer_aligned, st_smer_aligned_valid, st_mins, st_mins_valid);

    thread_dedup_v8<WINDOW, S_BITS>(st_mins, st_mins_valid, st_minimizers, st_minimizers_valid);

    thread_store_burst_v8<S_BITS>(st_minimizers, st_minimizers_valid, tab_hash, nMinizrs);
}
