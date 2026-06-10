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
void thr_min_v8(
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
        ap_uint<8> ocnt = 0;

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
