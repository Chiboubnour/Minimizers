#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
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
