#include <ap_int.h>
#include <hls_stream.h>

constexpr int SMER = 11;
constexpr int SMER_BITS = 2*SMER;

#if 0
//Bertrand
template<int SMER>
ap_uint<SMER> min2(ap_uint<SMER> a, ap_uint<SMER> b)
{
#pragma HLS INLINE
    return (a < b) ? a : b;
}

template<int SMER>
ap_uint<SMER> min8(ap_uint<SMER> tab[8])
{
#pragma HLS INLINE
    ap_uint<SMER> m1 = min2<SMER>(tab[0], tab[1]);
    ap_uint<SMER> m2 = min2<SMER>(tab[2], tab[3]);
    ap_uint<SMER> m3 = min2<SMER>(tab[4], tab[5]);
    ap_uint<SMER> m4 = min2<SMER>(tab[6], tab[7]);

    ap_uint<SMER> m5 = min2<SMER>(m1, m2);
    ap_uint<SMER> m6 = min2<SMER>(m3, m4);

    ap_uint<SMER> m7 = min2<SMER>(m5, m6);
    return m7;
}

template<int WINDOW,int SMER>
void thread_dedup(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> buf[WINDOW];
#pragma HLS ARRAY_PARTITION variable=buf complete

    ap_uint<SMER> last = ~ap_uint<SMER>(0);
    ap_uint<6> fill = 0;

MAIN_LOOP:
    while(true){
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<8*SMER> packed = in.read();
        ap_uint<8> valid = valid_i.read();

        if(valid == 0){
            out.write(0);
            valid_o.write(0);
            break;
        }

        ap_uint<SMER> mins[8];
        bool window_valid[8];
        ap_uint<8*SMER> pout = 0;
        ap_uint<4> ocnt = 0;

        ap_uint<6> fill_local = fill;

        // =========================
        // SHIFT + MIN
        // =========================
        for(int i=0;i<8;i++){
#pragma HLS UNROLL

            // lire la valeur
            ap_uint<SMER> v = packed.range((i+1)*SMER-1, i*SMER);
            ap_uint<SMER> v_masked = (i < valid) ? v : buf[0];

            // shift
            for(int j=WINDOW-1;j>0;j--){
#pragma HLS UNROLL
                buf[j] = buf[j-1];
            }
            buf[0] = v_masked;

            // update fill
            fill_local += (i < valid);
            bool win_valid = (fill_local >= WINDOW);

            // calcul du min sur la fenêtre
            ap_uint<SMER> m = buf[0];
            for(int j=1;j<WINDOW;j++){
#pragma HLS UNROLL
                m = (buf[j] < m) ? buf[j] : m;
            }

            mins[i] = m;
            window_valid[i] = (i < valid) && win_valid;
        }

        fill = fill_local;

        // =========================
        // DEDUP
        // =========================
        ap_uint<SMER> ref = last;
        ocnt = 0;


        for(int i=0;i<8;i++){
#pragma HLS UNROLL
            bool emit = window_valid[i] && (mins[i] != ref);
            ap_uint<SMER> val = mins[i];

            pout.range((ocnt+1)*SMER-1, ocnt*SMER) = emit ? val : pout.range((ocnt+1)*SMER-1, ocnt*SMER);
            ocnt += emit ? 1 : 0;


            ref = emit ? val : ref;
        }

        last = ref;

        if(ocnt > 0){
            out.write(pout);
            valid_o.write(ocnt);
        }
    }
}
#else
//Nour
template<int WINDOW,int SMER>
void thread_dedup(
    hls::stream<ap_uint<8*SMER>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER>>& out,
    hls::stream<ap_uint<8>>& valid_o
){
#pragma HLS INLINE off

    ap_uint<SMER> buf[WINDOW];
#pragma HLS ARRAY_PARTITION variable=buf complete

    ap_uint<SMER> last = ~ap_uint<SMER>(0);
    ap_uint<6> fill = 0;

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

        ap_uint<SMER> mins[8];

        bool window_valid[8];

        ap_uint<8*SMER> pout = 0;
        ap_uint<4> ocnt = 0;

        ap_uint<6> fill_local = fill;

        // =========================
        // SHIFT + MIN DIRECT
        // =========================
        for(int i=0;i<8;i++){

            bool active = (i < valid);

            ap_uint<SMER> v =
                packed.range((i+1)*SMER-1, i*SMER);

            ap_uint<SMER> v_masked =
                active ? v : buf[0];

            // SHIFT
            for(int j=WINDOW-1;j>0;j--){
                buf[j] = buf[j-1];
            }

            buf[0] = v_masked;

            // Mise à jour fill sans branchement
            fill_local += active;

            bool win_valid = (fill_local >= WINDOW);

            // Calcul MIN toujours exécuté
            ap_uint<SMER> m = buf[0];
            for(int j=1;j<WINDOW;j++){
                if(buf[j] < m)
                    m = buf[j];
            }

            mins[i] = m;
            window_valid[i] = active && win_valid;
        }

        fill = fill_local;

        // =========================
        // DEDUP
        // =========================
        ap_uint<SMER> ref = last;

        for(int i=0;i<8;i++){

            bool emit = window_valid[i] && (mins[i] != ref);

            ap_uint<SMER> m = mins[i];

            if(emit){
                pout.range((ocnt+1)*SMER-1, ocnt*SMER) = m;
                ocnt++;
            }

            if(window_valid[i]){
                ref = emit ? m : ref;
            }
        }

        last = ref;

        if(ocnt > 0){
            out.write(pout);
            valid_o.write(ocnt);
        }
    }
}
#endif

template void thread_dedup<4, 22>(
    hls::stream<ap_uint<8*22>>&,
    hls::stream<ap_uint<8>>&,
    hls::stream<ap_uint<8*22>>&,
    hls::stream<ap_uint<8>>&
);

