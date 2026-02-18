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

#if 0
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
#else
template <int S_BASE>
void thread_rolling_smer(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<2*S_BASE>>& smer_stream,
    hls::stream<bool>& smer_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;

    ap_uint<S_BITS> fwd = 0;
    ap_uint<S_BITS> rev = 0;
    ap_uint<16> base_cnt = 0;

    ap_uint<64> word = 0;
    ap_uint<8>  valid = 0;
    ap_uint<4>  idx = 8;   // force first read

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        // recharge word seulement quand nécessaire
        if(idx == 8){

            word  = base_stream.read();
            valid = base_valid.read();

            if(valid == 0){
                smer_valid.write(false);
                break;
            }

            idx = 0;
        }

        // lecture 1 base
        if(idx < valid){

            uint8_t nucl = (word >> (8*idx)) & 0xFF;
            ap_uint<2> c = (nucl >> 1) & 3;

            fwd = (fwd << 2) | c;
            rev = (rev >> 2) |
                  ((ap_uint<S_BITS>)(0x2 ^ c) << (S_BITS-2));

            base_cnt++;

            if(base_cnt >= S_BASE){

                ap_uint<S_BITS> canon =
                    (fwd < rev) ? fwd : rev;

                smer_stream.write(canon);
                smer_valid.write(true);
            }
        }

        idx++;
    }
}

template<int S_BITS>
void thread_hash_smer(
    hls::stream<ap_uint<S_BITS>>& smer_stream,
    hls::stream<bool>& smer_valid,

    hls::stream<ap_uint<S_BITS>>& hash_stream,
    hls::stream<bool>& hash_valid
){
#pragma HLS INLINE off

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        bool v = smer_valid.read();

        if(!v){
            hash_valid.write(false);
            break;
        }

        ap_uint<S_BITS> smer = smer_stream.read();

        ap_uint<S_BITS> h =
            hash_u64<S_BITS>((ap_uint<64>)smer);

        hash_stream.write(h);
        hash_valid.write(true);
    }
}
template<int S_BITS>
void thread_pack8(
    hls::stream<ap_uint<S_BITS>>& in,
    hls::stream<bool>& valid_i,

    hls::stream<ap_uint<8*S_BITS>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<8*S_BITS> buffer = 0;
    ap_uint<4> cnt = 0;

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        bool v = valid_i.read();

        if(!v){

            if(cnt > 0){
                out.write(buffer);
                valid_o.write(cnt);
            }

            out.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<S_BITS> h = in.read();

        buffer.range((cnt+1)*S_BITS-1, cnt*S_BITS) = h;
        cnt++;

        if(cnt == 8){
            out.write(buffer);
            valid_o.write(8);
            cnt = 0;
            buffer = 0;
        }
    }
}
#endif

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

    ap_uint<SMER> buf[WINDOW];
#pragma HLS ARRAY_PARTITION variable=buf complete

    ap_uint<SMER> last = ~ap_uint<SMER>(0);

    ap_uint<4> pos = 0;    
    ap_uint<5> fill = 0;

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

        ap_uint<8*SMER> pout = 0;
        ap_uint<4> ocnt = 0;

        for(int i=0;i<8;i++){

            if(i < valid){

                ap_uint<SMER> v =
                    packed.range((i+1)*SMER-1, i*SMER);

                buf[pos] = v;

                if(pos == WINDOW-1)
                    pos = 0;
                else
                    pos++;

                if(fill < WINDOW)
                    fill++;

                if(fill == WINDOW){


                    ap_uint<SMER> m0 = (buf[0]  < buf[1])  ? buf[0]  : buf[1];
                    ap_uint<SMER> m1 = (buf[2]  < buf[3])  ? buf[2]  : buf[3];
                    ap_uint<SMER> m2 = (buf[4]  < buf[5])  ? buf[4]  : buf[5];
                    ap_uint<SMER> m3 = (buf[6]  < buf[7])  ? buf[6]  : buf[7];
                    ap_uint<SMER> m4 = (buf[8]  < buf[9])  ? buf[8]  : buf[9];
                    ap_uint<SMER> m5 = (buf[10] < buf[11]) ? buf[10] : buf[11];
                    ap_uint<SMER> m6 = (buf[12] < buf[13]) ? buf[12] : buf[13];
                    ap_uint<SMER> m7 = (buf[14] < buf[15]) ? buf[14] : buf[15];

                    ap_uint<SMER> n0 = (m0 < m1) ? m0 : m1;
                    ap_uint<SMER> n1 = (m2 < m3) ? m2 : m3;
                    ap_uint<SMER> n2 = (m4 < m5) ? m4 : m5;
                    ap_uint<SMER> n3 = (m6 < m7) ? m6 : m7;

                    ap_uint<SMER> p0 = (n0 < n1) ? n0 : n1;
                    ap_uint<SMER> p1 = (n2 < n3) ? n2 : n3;

                    ap_uint<SMER> m  = (p0 < p1) ? p0 : p1;

                    if(m != last){
                        last = m;
                        pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                        ocnt++;
                    }
                }
            }
        }

        if(ocnt > 0){
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

    ap_uint<6> bcnt = 0;   

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

    hls::stream<ap_uint<S_BITS>> st_smer_roll;
    hls::stream<bool>            st_smer_roll_valid;

    hls::stream<ap_uint<S_BITS>> st_smer_hash;
    hls::stream<bool>            st_smer_hash_valid;

    hls::stream<ap_uint<8*S_BITS>> st_pack8;
    hls::stream<ap_uint<8>>        st_pack8_valid;

    hls::stream<ap_uint<8*S_BITS>> st_smer_aligned;
    hls::stream<ap_uint<8>>        st_smer_aligned_valid;

    hls::stream<ap_uint<8*S_BITS>> st_minimizers;
    hls::stream<ap_uint<8>>        st_minimizers_valid;

#pragma HLS STREAM variable=st_base_raw   depth=256
#pragma HLS STREAM variable=st_base_valid depth=256

#pragma HLS STREAM variable=st_base_aligned       depth=256
#pragma HLS STREAM variable=st_base_aligned_valid depth=256

#pragma HLS STREAM variable=st_smer_roll       depth=256
#pragma HLS STREAM variable=st_smer_roll_valid depth=256

#pragma HLS STREAM variable=st_smer_hash       depth=256
#pragma HLS STREAM variable=st_smer_hash_valid depth=256

#pragma HLS STREAM variable=st_pack8       depth=256
#pragma HLS STREAM variable=st_pack8_valid depth=256
   
#pragma HLS STREAM variable=st_smer_aligned       depth=256
#pragma HLS STREAM variable=st_smer_aligned_valid depth=256

#pragma HLS STREAM variable=st_minimizers       depth=256
#pragma HLS STREAM variable=st_minimizers_valid depth=256

    reader_hls( packed_sequence, n_bases,st_base_raw,st_base_valid );

    adapter_hls( st_base_raw, st_base_valid, st_base_aligned, st_base_aligned_valid);

    thread_rolling_smer<SMER>( st_base_aligned,st_base_aligned_valid, st_smer_roll,st_smer_roll_valid );

    thread_hash_smer<S_BITS>( st_smer_roll, st_smer_roll_valid, st_smer_hash, st_smer_hash_valid  );

    thread_pack8<S_BITS>( st_smer_hash, st_smer_hash_valid,st_pack8,st_pack8_valid );

    adapter_smer_hls<S_BITS>( st_pack8, st_pack8_valid, st_smer_aligned,st_smer_aligned_valid);

    thread_dedup_v8<WINDOW, S_BITS>(  st_smer_aligned, st_smer_aligned_valid, st_minimizers, st_minimizers_valid );

    thread_store_burst_v8<S_BITS>(  st_minimizers,st_minimizers_valid,tab_hash,nMinizrs );
}
