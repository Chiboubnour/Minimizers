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
void thr_dedup_v8(
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
        ap_uint<8> ocnt = 0;

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