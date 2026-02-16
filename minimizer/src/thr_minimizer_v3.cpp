#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include "functions.hpp"


constexpr int SMER   = 19;
constexpr int WINDOW = 16;

static void reader_hls(
    const uint64_t* packed_sequence,
    ap_uint<64> n_bases,
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid
) {
#pragma HLS INLINE off

    ap_uint<64> n_words = (n_bases + 7) >> 3;

    for (ap_uint<64> w = 0; w < n_words; w++) {

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
) {
#pragma HLS INLINE off

    ap_uint<128> buffer = 0;
    ap_uint<8>   valid_count = 0;


    while (true) {
#pragma HLS PIPELINE II=1

        ap_uint<64> in_word  = base_stream.read();
        ap_uint<8>  in_valid = base_valid.read();

        if (in_valid == 0) {
            // Écrire les s-mers restants dans le buffer
            if (valid_count > 0) {
                ap_uint<64> out_word = buffer.range(63,0);
                out_stream.write(out_word);
                out_valid.write(valid_count);
           
            }

            // Fin de flux
            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        // Ajouter le mot dans le buffer
        ap_uint<128> shifted = (ap_uint<128>)in_word << (valid_count * 8);
        buffer |= shifted;
        valid_count += in_valid;

        // Écrire tous les mots complets
        while (valid_count >= 8) {
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
) {
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * S_BASE;

    ap_uint<S_BITS> current_smer = 0;
    ap_uint<S_BITS> cur_inv_smer = 0;
    int base_count = 0;

    ap_uint<8*S_BITS> buffer = 0;
    ap_uint<8> count_smers = 0;


    while (true) {
#pragma HLS PIPELINE II=1
        ap_uint<64> word  = base_stream.read();
        ap_uint<8>  valid = base_valid.read();

        if (valid == 0) {
            // Écrire uniquement les SMER valides restants
            if (count_smers > 0) {
                out_stream.write(buffer);
                out_valid.write(count_smers);
            }
            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        for (int i = 0; i < valid; i++) {
#pragma HLS UNROLL factor=1
            uint8_t nucl = (word >> (8*i)) & 0xFF;
            uint8_t c_nucl = (nucl >> 1) & 0x03;

            // Rolling forward
            current_smer <<= 2;
            current_smer(1,0) = c_nucl;
            current_smer &= ((ap_uint<S_BITS>) - 1);

            // Rolling reverse
            cur_inv_smer >>= 2;
            cur_inv_smer(S_BITS-1, S_BITS-2) = 0x2 ^ c_nucl;

            base_count++;

            if (base_count >= S_BASE) {
                ap_uint<S_BITS> vmin  = min_v1<S_BITS>(current_smer, cur_inv_smer);
                ap_uint<S_BITS> vhash = hash_u64<S_BITS>(vmin);

                buffer.range(count_smers*S_BITS + S_BITS-1,
                             count_smers*S_BITS) = vhash;
                count_smers++;

                if (count_smers == 8) {
                    out_stream.write(buffer);
                    out_valid.write(8);
                    buffer = 0;
                    count_smers = 0;
                }
            }
        }
    }

    if (count_smers > 0) {
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

    ap_uint<8*S_BITS> buffer = 0;
    ap_uint<8> valid_count = 0;

    unsigned iter = 0;
    unsigned total_out = 0;

    while(true){
#pragma HLS PIPELINE II=1

        iter++;

        auto in_word  = base_stream.read();
        auto in_valid = base_valid.read();

        // ===== EOS =====
        if(in_valid == 0){

            if(valid_count > 0){

                out_stream.write(buffer);
                out_valid.write(valid_count);

                total_out += valid_count;
            }

            out_stream.write(0);
            out_valid.write(0);
            break;
        }

        // ===== CONCAT =====
        buffer.range((valid_count+in_valid)*S_BITS-1,
                     valid_count*S_BITS) =
            in_word.range(in_valid*S_BITS-1, 0);

        valid_count += in_valid;

        // ===== PACKETS DE 8 =====
        while(valid_count >= 8){

            ap_uint<8*S_BITS> out_word =
                buffer.range(8*S_BITS-1,0);

            out_stream.write(out_word);
            out_valid.write(8);

            total_out += 8;

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
    int pos = 0;
    int fill = 0; 

    unsigned iter = 0;
    unsigned total_out = 0;

    while(true){
#pragma HLS PIPELINE II=1

        iter++;

        auto packed = in.read();
        auto valid  = valid_i.read();

        if(valid == 0){
            out.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<8*SMER> pout = 0;
        ap_uint<8> ocnt = 0;

        for(int i=0;i<valid;i++){
#pragma HLS UNROLL
            ap_uint<SMER> v = packed.range((i+1)*SMER-1, i*SMER);

            buf[pos] = v;
            pos = (pos+1)%WINDOW;

            if(fill < WINDOW) {
                fill++; // remplir la fenêtre
            }

            // Calculer le minimizer dès que fill >= WINDOW
            if(fill >= WINDOW){
                ap_uint<SMER> m = buf[0];
                for(int j=1;j<WINDOW;j++){
#pragma HLS UNROLL
                    if(buf[j] < m) m = buf[j];
                }

                if(m != last){
                    last = m;
                    pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                    ocnt++;
                    total_out++;
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

    while (true) {
#pragma HLS PIPELINE II=1
        auto packed = in.read();
        auto valid  = valid_i.read();

        if (valid == 0) break;

        for (int i = 0; i < valid; i++) {
#pragma HLS UNROLL
            ap_uint<SMER> v = packed.range((i+1)*SMER-1, i*SMER);
            tab_hash[cnt++] = v;   
        }
    }

    *nElem = cnt; 
}



void minimizer(
    const uint64_t* seq,
    ap_uint<512>* tab_hash,
    ap_uint<64>* nMin,
    ap_uint<64> n_bases
){
#pragma HLS DATAFLOW

    constexpr int S_BITS = 2 * SMER;

    hls::stream<ap_uint<64>> s0,s1;
    hls::stream<ap_uint<8>>  v0,v1;

    hls::stream<ap_uint<8*S_BITS>> s2,s3;
    hls::stream<ap_uint<8>> v2,v3;

    hls::stream<ap_uint<8*S_BITS>> s4;
    hls::stream<ap_uint<8>> v4;

    reader_hls(seq, n_bases, s0, v0);
    adapter_hls(s0, v0, s1, v1);

    smer_thread_hls<SMER>(s1, v1, s2, v2);

    adapter_smer_hls<S_BITS>(s2, v2, s3, v3);

    thread_dedup_v8<WINDOW, S_BITS>(s3, v3, s4, v4);

    thread_store_burst_v8<S_BITS>(s4, v4, reinterpret_cast<ap_uint<64>*>(tab_hash), nMin);
}


extern "C"
void minimizer_v3(
    const uint64_t* seq,
    int s,int w,int n,
    uint64_t* tab_hash,
    uint64_t* nMin
){
    if(s != 19 || w != 16){
        *nMin = 0;
        return;
    }

    ap_uint<512>* th = (ap_uint<512>*)tab_hash;
    ap_uint<64> nm = 0;

    minimizer(seq, th, &nm, n);

    *nMin = nm;
}
