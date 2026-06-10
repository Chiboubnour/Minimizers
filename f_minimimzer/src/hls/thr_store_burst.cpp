#pragma once
#include "./header.hpp"
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
void thr_store_burst(
    hls::stream< ap_uint<PAR_FACTOR * SMER_SIZE> >& minz_i,
    hls::stream< ap_uint<PAR_FACTOR>             >& minz_v,
    ap_uint<64>* tab_hash,
    ap_uint<64>* nElem
){
#pragma HLS INLINE off

    const ap_uint<SMER_SIZE> zero_s = 0;
    const ap_uint<        1> zero_v = 0;
    //
    //
    ap_uint<64> cnt = 0;
    //
    //
    ap_uint<PAR_FACTOR * SMER_SIZE> packed_v;
    ap_uint<PAR_FACTOR            > valid_v;
    //
    //    
    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        //
        //
        packed_v = minz_i.read();
        valid_v  = minz_v.read();
        //
        //
        if (valid_v == 0){
            break;
        }
        //
        //
        for(int i = 0; i < PAR_FACTOR; i++)
        {
            //
            //
            const ap_uint<SMER_SIZE> final_m = packed_v.range(SMER_SIZE-1, 0);
            if(valid_v[0] == 1) {
                tab_hash[cnt++] = final_m;
            }
            //
            //
            const ap_uint<(PAR_FACTOR-1) * SMER_SIZE> up_part = packed_v.range(PAR_FACTOR * SMER - 1, SMER);
            packed_v = (zero_s, up_v);
            //
            //
        }
    }
    //
    //
    *nElem = cnt;
    //
    //
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