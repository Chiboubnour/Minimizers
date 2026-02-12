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

    for (ap_uint<64> w=0; w<n_words; w++) {

        ap_uint<64> word = packed_sequence[w];
        ap_uint<8> valid;
        if(w == n_words-1 && (n_bases & 7)) {
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

#if 0
static void adapter_bases_v8(
    hls::stream<ap_uint<64>>& base_i,
    hls::stream<ap_uint<8>>&  valid_i,
    hls::stream<ap_uint<64>>& base_o,
    hls::stream<ap_uint<8>>&  valid_o
){
#pragma HLS INLINE off

    while(true){
#pragma HLS PIPELINE II=1

        auto b = base_i.read();
        auto v = valid_i.read();

        if(v==0){
            base_o.write(0);
            valid_o.write(0);
            break;
        }

        base_o.write(b);
        valid_o.write(v);
    }
}
#else
static void adapter_bases_v8(
    hls::stream<ap_uint<64>>& base_i,
    hls::stream<ap_uint<8>>&  valid_i,
    hls::stream<ap_uint<64>>& base_o,
    hls::stream<ap_uint<8>>&  valid_o
){
#pragma HLS INLINE off

    while(true){
#pragma HLS PIPELINE II=1

        const ap_uint<64> b = base_i.read();
        const ap_uint<8> v = valid_i.read();

        printf("[adapter_bases_v8] Read word   = 0x%016llx, valid = 0x%02x\n",
               (unsigned long long)b, (unsigned int)v);

        if(v != 0xFF){
            base_o.write(0);
            valid_o.write(0);
            printf("[adapter_bases_v8] End of stream, 0 written\n");
            break;
        }

        base_o.write(b);
        valid_o.write(v);

        printf("[adapter_bases_v8] Written word = 0x%016llx, valid = 0x%02x\n",
               (unsigned long long)b, (unsigned int)v);
    }
}
#endif


#if 0
template<int SMER>
void thread_smer_generator_v2(
    hls::stream<ap_uint<64>>& base_i,
    hls::stream<ap_uint<8>>&  valid_i,
    hls::stream<ap_uint<8*SMER>>& smer_o,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    constexpr int INNER = 2*SMER;

    ap_uint<INNER> fwd = 0;
    ap_uint<INNER> rev  = 0;
    ap_uint<6> init     = 0;
    bool ready          = false;

    while(true){
#pragma HLS PIPELINE II=1

        auto bases = base_i.read();
        auto valid = valid_i.read();

        if(valid==0){
            smer_o.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<8*SMER> pack = 0;
        ap_uint<8> cnt = 0;

        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            if(!valid[i]) continue;

            ap_uint<8> n = bases.range(8*i+7,8*i);
            if(n==0) continue;

            ap_uint<2> c = (n>>1)&3;

            fwd <<= 2;
            rev >>= 2;
            fwd(1,0) = c;
            rev(INNER-1,INNER-2) = 0x2 ^ c;

            if(!ready){
                init++;
                if(init >= SMER-1) ready = true;
            }

            if(ready){
                auto vmin = (fwd<rev)?fwd:rev;
                ap_uint<SMER> hv = hash_u64(vmin).range(SMER-1,0);

                pack.range((cnt+1)*SMER-1,cnt*SMER) = hv;
                cnt++;
            }
        }

        if(cnt > 0){
            smer_o.write(pack);
            valid_o.write(cnt);
        } else {
            printf("[SMER_GEN] NO smer generated for this word\n");
        }
    }
}

#else
template<int SMER>
void thread_smer_generator_v2(
    hls::stream<ap_uint<64>>& base_i,
    hls::stream<ap_uint<8>>&  valid_i,
    hls::stream<ap_uint<8*SMER>>& smer_o,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    constexpr int INNER = 2*SMER; 
    ap_uint<INNER> fwd = 0;      
    ap_uint<INNER> rev = 0;      
    ap_uint<6> filled = 0;          

    while(true){
#pragma HLS PIPELINE II=1

        auto bases = base_i.read();
        auto valid = valid_i.read();

        printf("[SMER_GEN] Read word: 0x%016llx valid=%u\n",
               (unsigned long long)bases, (unsigned)valid);

        if(valid == 0){               // fin du stream
            printf("[SMER_GEN] EOS received\n");
            smer_o.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<8*SMER> pack = 0;     // buffer pour stocker les S-mers générés
        ap_uint<8> cnt = 0;

        for(int i=0; i<8; i++){
#pragma HLS UNROLL

            const ap_uint<8> n = bases.range(8*i+7, 8*i);
            const ap_uint<2> c = (n >> 1) & 0x3;
            if(n == 0){
                printf("%s:%d : thread_smer_generator_v2[%d] n = %d and c = %d\n", __FILE__, __LINE__, i, n.to_int(), c.to_int());
                continue;
            }      

            // mise à jour des S-mers forward et reverse
            fwd <<= 2;
            fwd(1,0) = c;
            rev >>= 2;
            rev(INNER-1, INNER-2) = 0x2 ^ c;  // complément inversé

            filled++;                   // base accumulée

            printf("  [SMER_GEN] base[%d]=%u encoded=%u fwd=0x%llx rev=0x%llx filled=%u\n",
                   i, (unsigned)n, (unsigned)c,
                   (unsigned long long)fwd,
                   (unsigned long long)rev,
                   (unsigned)filled);

            if(filled >= SMER){
                ap_uint<INNER> vmin = (fwd < rev) ? fwd : rev;
                ap_uint<SMER> hv = hash_u64<SMER>(vmin);  

                pack.range(cnt*SMER + SMER -1, cnt*SMER) = hv;
                cnt++;

                printf("    [SMER_GEN] Generated S-mer %u: vmin=0x%llx hv=0x%llx\n",
                       (unsigned)cnt,
                       (unsigned long long)vmin,
                       (unsigned long long)hv);
            }
        }

        if(cnt > 0){
            printf("[SMER_GEN] Writing %u S-mers\n", (unsigned)cnt);
            smer_o.write(pack);
            valid_o.write(cnt);
        } else {
            printf("[SMER_GEN] No S-mer generated for this word\n");
        }
    }
}

#endif

#if 0
template<int SMER>
void adapter_smers_v8(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> fifo[8];
#pragma HLS ARRAY_PARTITION variable=fifo complete

    ap_uint<4> fill=0;

    while(true){
#pragma HLS PIPELINE II=1

        auto packed = in.read();
        auto valid  = valid_i.read();

        if(valid==0){

            if(fill>0){
                ap_uint<8*SMER> w=0;
                for(int i=0;i<fill;i++){
#pragma HLS UNROLL
                    w.range((i+1)*SMER-1,i*SMER)=fifo[i];
                }
                out.write(w);
                valid_o.write(fill);
            }

            out.write(0);
            valid_o.write(0);
            break;
        }

        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            if(i<valid){
                fifo[fill++] = packed.range((i+1)*SMER-1,i*SMER);

                if(fill==8){
                    ap_uint<8*SMER> w=0;
                    for(int j=0;j<8;j++){
#pragma HLS UNROLL
                        w.range((j+1)*SMER-1,j*SMER)=fifo[j];
                    }

                    out.write(w);
                    valid_o.write(8);
                    fill=0;
                }
            }
        }
    }
}

#else
template<int SMER>
void adapter_smers_v8(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> fifo[8];
#pragma HLS ARRAY_PARTITION variable=fifo complete

    ap_uint<4> fill=0;

    while(true){
#pragma HLS PIPELINE II=1

        auto packed = in.read();
        auto valid  = valid_i.read();

        printf("[adapter_smers_v8] Read packed, valid=%u, fill=%u\n",
               (unsigned)valid, (unsigned)fill);

        if(valid==0){

            if(fill>0){
                ap_uint<8*SMER> w=0;

                for(int i=0;i<fill;i++){
#pragma HLS UNROLL
                    w.range((i+1)*SMER-1,i*SMER)=fifo[i];
                }

                out.write(w);
                valid_o.write(fill);

                printf("[EOS flush] fill=%u\n",(unsigned)fill);
            }

            out.write(0);
            valid_o.write(0);
            break;
        }

        for(int i=0;i<valid;i++){
#pragma HLS UNROLL

            fifo[fill++] = packed.range((i+1)*SMER-1,i*SMER);

            printf("  [fifo add] fill=%u\n",(unsigned)fill);

            if(fill==8){

                ap_uint<8*SMER> w=0;

                for(int j=0;j<8;j++){
#pragma HLS UNROLL
                    w.range((j+1)*SMER-1,j*SMER)=fifo[j];
                }

                out.write(w);
                valid_o.write(8);

                printf("[FULL flush]\n");

                fill=0;
            }
        }
    }
}


#endif

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

    // Initialisation buffer avec valeur max (sentinelle)
    for(int i=0;i<WINDOW;i++){
#pragma HLS UNROLL
        buf[i] = ~ap_uint<SMER>(0);
    }

    ap_uint<SMER> last = ~ap_uint<SMER>(0);
    int pos = 0;

    while(true){
#pragma HLS PIPELINE II=1

        auto packed = in.read();
        auto valid  = valid_i.read();

        // Fin de stream
        if(valid == 0){
            out.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<8*SMER> pout = 0;
        ap_uint<8> ocnt = 0;

        // ----------- traitement paquet -----------
        for(int i=0;i<valid;i++){
#pragma HLS UNROLL

            ap_uint<SMER> v =
                packed.range((i+1)*SMER-1, i*SMER);

            // Mise à jour fenêtre circulaire
            buf[pos] = v;

            // Recherche minimizer
            ap_uint<SMER> m = buf[0];
            for(int j=1;j<WINDOW;j++){
#pragma HLS UNROLL
                if(buf[j] < m)
                    m = buf[j];
            }

            pos++;
            if(pos == WINDOW)
                pos = 0;

            // Deduplication
            if(m != last){
                last = m;

                pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                ocnt++;
            }
        }

        // Écriture sortie si données
        if(ocnt > 0){
            out.write(pout);
            valid_o.write(ocnt);
        }
    }
}


template<int SMER>
void thread_store_burst_v8(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    ap_uint<512>* tab_hash,
    ap_uint<64>* nElem
){
#pragma HLS INLINE off

    ap_uint<64> cnt=0;
    ap_uint<64> mem=0;

    while(true){
#pragma HLS PIPELINE II=1

        auto packed = in.read();
        auto valid  = valid_i.read();

        if(valid==0) break;

        ap_uint<512> w=0;

        for(int i=0;i<valid;i++){
#pragma HLS UNROLL

            ap_uint<64> v = 0;

            v.range(SMER-1,0) =
                packed.range((i+1)*SMER-1,i*SMER);

            w.range((i+1)*64-1,i*64) = v;
        }

        tab_hash[mem++] = w;
        cnt += valid;
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

    hls::stream<ap_uint<64>> s0,s1;
    hls::stream<ap_uint<8>>  v0,v1;

    hls::stream<ap_uint<8*SMER>> s2,s3;
    hls::stream<ap_uint<8>> v2,v3;

    hls::stream<ap_uint<8*SMER>> s4;
    hls::stream<ap_uint<8>> v4;

#pragma HLS STREAM variable=s0 depth=256
#pragma HLS STREAM variable=s1 depth=256
#pragma HLS STREAM variable=s2 depth=256
#pragma HLS STREAM variable=s3 depth=256
#pragma HLS STREAM variable=s4 depth=256

#pragma HLS STREAM variable=v0 depth=256
#pragma HLS STREAM variable=v1 depth=256
#pragma HLS STREAM variable=v2 depth=256
#pragma HLS STREAM variable=v3 depth=256
#pragma HLS STREAM variable=v4 depth=256

    reader_hls(seq,n_bases,s0,v0);
    adapter_bases_v8(s0,v0,s1,v1);
    thread_smer_generator_v2<SMER>(s1,v1,s2,v2);
    adapter_smers_v8<SMER>(s2,v2,s3,v3);
    thread_dedup_v8<WINDOW,SMER>(s3,v3,s4,v4);
    thread_store_burst_v8<SMER>(s4,v4,tab_hash,nMin);
}

extern "C"
void minimizer_v3(
    const uint64_t* seq,
    int s,int w,int n,
    uint64_t* tab_hash,
    uint64_t* nMin
){
    if(s!=19 || w!=16){
        *nMin=0;
        return;
    }

    ap_uint<512>* th = (ap_uint<512>*)tab_hash;
    ap_uint<64> nm=0;

    minimizer(seq,th,&nm,n);

    *nMin=nm;
}
