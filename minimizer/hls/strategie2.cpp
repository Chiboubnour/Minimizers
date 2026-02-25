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

    // rolling smer
    ap_uint<S_BITS> current_smer = 0;
    ap_uint<S_BITS> cur_inv_smer  = 0;
    int base_count = 0;

    // buffer pour 8 smers avant écriture
    ap_uint<8*S_BITS> buffer = 0;
    ap_uint<8> count_smers = 0;

    bool init_done = false;

    // ----- phase principale -----
    while (true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if(valid == 0) break;

        // ---- décodage des 8 bases ----
        for(int i=0;i<8;i++){
            if(i >= valid) break;

            uint8_t nucl = (word >> (8*i)) & 0xFF;
            ap_uint<2> c_nucl = (nucl >> 1) & 0x3;

            // ----- phase d'initialisation -----
            if(!init_done){
                current_smer <<= 2;
                current_smer(1,0) = c_nucl;

                cur_inv_smer >>= 2;
                cur_inv_smer(S_BITS-1,S_BITS-2) = 0x2 ^ c_nucl;

                base_count++;
                if(base_count == S_BASE){
                    init_done = true;
                    // premier smer généré
                    ap_uint<S_BITS> vmin  =
                        (current_smer < cur_inv_smer) ? current_smer : cur_inv_smer;
                    ap_uint<S_BITS> vhash = hash_u64<S_BITS>((ap_uint<64>)vmin);

                    buffer.range((count_smers+1)*S_BITS-1, count_smers*S_BITS) = vhash;
                    count_smers++;
                }
                continue;
            }

            // ----- phase rolling -----
            current_smer <<= 2;
            current_smer(1,0) = c_nucl;

            cur_inv_smer >>= 2;
            cur_inv_smer(S_BITS-1,S_BITS-2) = 0x2 ^ c_nucl;

            // rolling smer
            ap_uint<S_BITS> vmin  =
                (current_smer < cur_inv_smer) ? current_smer : cur_inv_smer;
            ap_uint<S_BITS> vhash = hash_u64<S_BITS>((ap_uint<64>)vmin);

            buffer.range((count_smers+1)*S_BITS-1, count_smers*S_BITS) = vhash;
            count_smers++;

            // écriture tous les 8 smers
            if(count_smers == 8){
                out_stream.write(buffer);
                out_valid.write(8);
                buffer = 0;
                count_smers = 0;
            }
        }
    }

    if(count_smers > 0){
        out_stream.write(buffer);
        out_valid.write(count_smers);
    }

    out_stream.write(0);
    out_valid.write(0);
}
#else
template<int S_BASE>
void smer_thread_hls(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<8*2*S_BASE>>& out_stream,
    hls::stream<ap_uint<8>>&  out_valid
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2*S_BASE;
    const ap_uint<S_BITS> MASK = ap_uint<S_BITS>(-1);

    ap_uint<S_BITS> fwd_state = 0;
    ap_uint<S_BITS> rev_state = 0;
    ap_uint<32> base_cnt = 0;

    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if(valid == 0){
            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        ap_uint<8*S_BITS> packed_out = 0;
        ap_uint<4> out_count = 0;

        ap_uint<S_BITS> fwd_local = fwd_state;
        ap_uint<S_BITS> rev_local = rev_state;

        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            if(i < valid){

                uint8_t nucl = (word >> (8*i)) & 0xFF;
                ap_uint<2> c = (nucl >> 1) & 3;
                ap_uint<2> comp = 0x2 ^ c;

                fwd_local <<= 2;
                fwd_local(1,0) = c;

                rev_local >>= 2;
                rev_local(S_BITS-1,S_BITS-2) = comp;

                base_cnt++;

                if(base_cnt >= S_BASE){

                    ap_uint<S_BITS> canon = (fwd_local < rev_local) ? fwd_local : rev_local;

                    ap_uint<S_BITS> h = hash_u64<S_BITS>((ap_uint<64>)canon);

                    packed_out.range( (out_count+1)*S_BITS-1,  out_count*S_BITS ) = h;

                    out_count++;
                }
            }
        }

        fwd_state = fwd_local & MASK;
        rev_state = rev_local;

        if(out_count > 0){
            out_stream.write(packed_out);
            out_valid.write(out_count);
        }
    }
}
#endif
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
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
#if 1
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
    ap_uint<4> pos  = 0;
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
#pragma HLS LOOP_FLATTEN off

            if(i >= valid) break;

#pragma HLS PIPELINE II=1

            ap_uint<SMER> v = packed.range((i+1)*SMER-1, i*SMER);

            // insertion circulaire
            buf[pos] = v;
            pos = (pos + 1) % WINDOW;

            if(fill < WINDOW) fill++;

            // fenêtre valide ?
            if(fill == WINDOW){
                ap_uint<SMER> m = buf[0];
                for(int j=1;j<WINDOW;j++){
#pragma HLS UNROLL
                    if(buf[j] < m) m = buf[j];
                }

                // sortie si nouveau minimum
                if(m != last){
                    last = m;
                    pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                    ocnt++;
                }
            }
        }

        if(ocnt > 0){
            out.write(pout);
            valid_o.write(ocnt);
        }
    }
}
#else
template<int WINDOW, int SMER>
void packet_window_min_v8(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> buf[WINDOW];
#pragma HLS ARRAY_PARTITION variable=buf complete

    ap_uint<4> pos  = 0;
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
#pragma HLS UNROLL

            if(i >= valid) break;

            ap_uint<SMER> v =
                packed.range((i+1)*SMER-1, i*SMER);

            // insertion circulaire
            buf[pos] = v;
            pos = (pos + 1) % WINDOW;

            if(fill < WINDOW) fill++;

            if(fill == WINDOW){

                // calcul min exact sur 16
                ap_uint<SMER> m = buf[0];
                for(int j=1;j<WINDOW;j++){
#pragma HLS UNROLL
                    if(buf[j] < m)
                        m = buf[j];
                }

                pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                ocnt++;
            }
        }

        if(ocnt > 0){
            out.write(pout);
            valid_o.write(ocnt);
        }
    }
}

template<int SMER>
void packet_dedup_v8(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& out,
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

        ap_uint<8*SMER> pout = 0;
        ap_uint<4> ocnt = 0;

        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            if(i >= valid) break;

            ap_uint<SMER> v =
                packed.range((i+1)*SMER-1, i*SMER);

            if(v != last){
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
#endif
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


    reader_hls(  packed_sequence,n_bases,st_base_raw,st_base_valid );

    adapter_hls(st_base_raw, st_base_valid, st_base_aligned,  st_base_aligned_valid );

    smer_thread_hls<SMER>( st_base_aligned, st_base_aligned_valid, st_smer_packed,st_smer_packed_valid);

    adapter_smer_hls<S_BITS>( st_smer_packed, st_smer_packed_valid,st_smer_aligned,st_smer_aligned_valid);
#if 1
    thread_dedup_v8<WINDOW, S_BITS>( st_smer_aligned,st_smer_aligned_valid,st_minimizers,st_minimizers_valid);
#else

    hls::stream<ap_uint<8*S_BITS>> st_window_min;
    hls::stream<ap_uint<8>>        st_window_min_valid;

#pragma HLS STREAM variable=st_window_min       depth=256
#pragma HLS STREAM variable=st_window_min_valid depth=256

    packet_window_min_v8<WINDOW, S_BITS>( st_smer_aligned,st_smer_aligned_valid,st_window_min, st_window_min_valid);
    packet_dedup_v8<S_BITS>(st_window_min,st_window_min_valid,st_minimizers, st_minimizers_valid);

#endif
    thread_store_burst_v8<S_BITS>(st_minimizers,st_minimizers_valid,tab_hash,nMinizrs);
}
