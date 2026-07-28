#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include "minimizer.hpp"

//
//
//
//
template<int PAR_FACTOR>
void thr_axis_unpack(
    hls::stream<ap_uint<AXIS_BITS>>&       in_axis,
    hls::stream<ap_uint<64 * PAR_FACTOR>>& out_stream,
    hls::stream<ap_uint<8>>&               out_count,
    hls::stream<bool>&                     out_last
){
#pragma HLS INLINE off
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<AXIS_BITS> packed = in_axis.read();
        const ap_uint<64 * PAR_FACTOR> word  = packed.range(64 * PAR_FACTOR - 1, 0);
        const ap_uint<8>               count = packed.range(AXIS_BITS - 1, 64 * PAR_FACTOR);
        const bool                     last  = (count == 0);

        out_stream.write(word);
        out_count.write(count);
        out_last.write(last);
        if (last) break;
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
template<int PAR_FACTOR>
void thr_flatten(
    hls::stream<ap_uint<64 * PAR_FACTOR>>& elem_stream_i,
    hls::stream<ap_uint<8>>&               elem_count_i,
    hls::stream<bool>&                     elem_last_i,
    hls::stream<ap_uint<GROUP_BITS>>&      out_v,
    hls::stream<bool>&                     out_e,
    hls::stream<ap_uint<64>>&              total_o)
{
#pragma HLS INLINE off

    ap_uint<GROUP_BITS> group_buf = 0;
    ap_uint<8>          group_pos = 0;   // 0..GROUP_W-1
    ap_uint<64>         total     = 0;


    ap_uint<64 * PAR_FACTOR> word    = elem_stream_i.read();
    ap_uint<8>               remain  = elem_count_i.read();
    bool                     is_end0 = elem_last_i.read();
    ap_uint<8>               pos     = 0;

    if (is_end0) {
        out_v.write(0);
        out_e.write(false);   // sentinelle de fin
        total_o.write(total);
        return;
    }

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<64> elem = word.range(64 * (pos + 1) - 1, 64 * pos);
        pos++;
        total++;

        group_buf.range(64 * (group_pos + 1) - 1, 64 * group_pos) = elem;
        group_pos++;

        if (group_pos == GROUP_W) {
            out_v.write(group_buf);
            out_e.write(true);
            group_buf = 0;
            group_pos = 0;
        }

        if (remain == 1) {

            word = elem_stream_i.read();
            const ap_uint<8> count  = elem_count_i.read();
            const bool       is_end = elem_last_i.read();

            if (is_end) {
                if (group_pos != 0) {
                    out_v.write(group_buf);
                    out_e.write(true);
                }
                out_v.write(0);
                out_e.write(false);   // sentinelle de fin
                total_o.write(total);
                break;
            }
            remain = count;
            pos    = 0;
        } else {
            remain--;
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
void thr_write_burst(
    hls::stream<ap_uint<GROUP_BITS>>& inv,
    hls::stream<bool>&                ine,
    hls::stream<ap_uint<64>>&          total_i,
    ap_uint<GROUP_BITS>*               dst,
    ap_uint<64>*                       nMinizrs,
    ap_uint<64>                        n_bases)
{
#pragma HLS INLINE off

    bool                lastE = false;
    ap_uint<GROUP_BITS> x     = 0;
    bool                e     = false;

    const ap_uint<64> n_groups_max = (n_bases + GROUP_W - 1) / GROUP_W;
    for (ap_uint<64> i = 0; i < n_groups_max; i += BURST_LEN)
    {
#pragma HLS LOOP_TRIPCOUNT min=64 max=64 avg=64
#pragma HLS DEPENDENCE variable=dst class=array inter false
        for (int j = 0; j < BURST_LEN; j++)
        {
#pragma HLS PIPELINE II=1
            if (!lastE) {
                x = inv.read();
                e = ine.read();
                if (!e) lastE = true;
            }
            dst[i + j] = x;   // ecriture inconditionnelle (hors du if) --
                               // necessaire pour l'inference de rafale m_axi.
        }
        if (lastE) break;
    }


    *nMinizrs = total_i.read();
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
void minimizer_store(
    hls::stream<ap_uint<AXIS_BITS>>& elem_axis_i,
    ap_uint<GROUP_BITS>*              tab_hash,
    ap_uint<64>*                      nMinizrs,
    ap_uint<64>                       n_bases
){
#pragma HLS INTERFACE axis      port=elem_axis_i
#pragma HLS INTERFACE m_axi     port=tab_hash offset=slave bundle=gmem1 num_write_outstanding=32 depth=COSIM_MAX_HASHGROUP max_write_burst_length=128
#pragma HLS INTERFACE s_axilite port=nMinizrs
#pragma HLS INTERFACE s_axilite port=n_bases
#pragma HLS INTERFACE s_axilite port=return

#pragma HLS DATAFLOW
    hls::stream<ap_uint<64 * PAR_FACTOR>> s_elem;
    hls::stream<ap_uint<8>>               s_elem_count;
    hls::stream<bool>                     s_elem_last;
    hls::stream<ap_uint<GROUP_BITS>>      s_flat_v;
    hls::stream<bool>                     s_flat_e;
    hls::stream<ap_uint<64>>              s_total;

#pragma HLS STREAM variable=s_elem       depth=64
#pragma HLS STREAM variable=s_elem_count depth=64
#pragma HLS STREAM variable=s_elem_last  depth=64
#pragma HLS STREAM variable=s_flat_v     depth=64
#pragma HLS STREAM variable=s_flat_e     depth=64
#pragma HLS STREAM variable=s_total      depth=2

    thr_axis_unpack<PAR_FACTOR>(elem_axis_i, s_elem, s_elem_count, s_elem_last);
    thr_flatten<PAR_FACTOR>(s_elem, s_elem_count, s_elem_last, s_flat_v, s_flat_e, s_total);
    thr_write_burst(s_flat_v, s_flat_e, s_total, tab_hash, nMinizrs, n_bases);
}
