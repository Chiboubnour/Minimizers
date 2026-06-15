#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m6 : aligneur de fenetre (template sur PAR_FACTOR), header-only.
//
//  Meme operation que m3, mais warm-up de (WINDOW_SIZE-1) elements, sur des hashes
//  de 2*SMER_SIZE bits. Re-decoupe le flux dense en first_rounds paquets pleins +
//  1 paquet de n_last_round, puis des paquets pleins.
//
//      first_rounds = (WINDOW_SIZE-1) / PAR_FACTOR
//      n_last_round = (WINDOW_SIZE-1) % PAR_FACTOR
//      buff_values  = PAR_FACTOR - n_last_round
//
template<int PAR_FACTOR>
void thr_adapter_smer(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS       = 2 * SMER_SIZE;                      // un hash
    constexpr int first_rounds = (WINDOW_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (WINDOW_SIZE - 1) % PAR_FACTOR;
    constexpr int buff_values  = PAR_FACTOR - n_last_round;
    const ap_uint<PAR_FACTOR> last_valid = ((ap_uint<PAR_FACTOR>)1 << n_last_round) - 1;

    ap_uint<buff_values * S_BITS> buffer_d = 0;
    ap_uint<buff_values>          buffer_v = 0;

    ap_uint<S_BITS * PAR_FACTOR> in_word;
    ap_uint<         PAR_FACTOR> in_valid;

    // first_rounds paquets pleins verbatim
    for (int r = 0; r < first_rounds; r++) {
        in_word  = smer_stream_i.read();
        in_valid = smer_valid_i.read();
        smer_stream_o.write(in_word);
        smer_valid_o.write (in_valid);
    }
    // paquet de warm-up partiel
    in_word  = smer_stream_i.read();
    in_valid = smer_valid_i.read();
    smer_stream_o.write(in_word);
    smer_valid_o.write (last_valid);
    buffer_d = in_word.range (S_BITS * PAR_FACTOR - 1, S_BITS * n_last_round);
    buffer_v = in_valid.range(         PAR_FACTOR - 1,          n_last_round);

    // regime permanent
    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        in_word  = smer_stream_i.read();
        in_valid = smer_valid_i.read();

        const ap_uint<buff_values  * S_BITS> up_part_d = in_word.range (S_BITS * PAR_FACTOR - 1, S_BITS * n_last_round);
        const ap_uint<buff_values>           up_part_v = in_valid.range(         PAR_FACTOR - 1,          n_last_round);
        const ap_uint<n_last_round * S_BITS> dw_part_d = in_word.range (S_BITS * n_last_round - 1, 0);
        const ap_uint<n_last_round>          dw_part_v = in_valid.range(         n_last_round - 1, 0);

        const ap_uint<S_BITS * PAR_FACTOR> to_send_d = (dw_part_d, buffer_d);
        const ap_uint<         PAR_FACTOR> to_send_v = (dw_part_v, buffer_v);

        if (to_send_v != 0) {
            smer_stream_o.write(to_send_d);
            smer_valid_o.write (to_send_v);
        }

        buffer_d = up_part_d;
        buffer_v = up_part_v;

        if (in_valid == 0) {
            break;
        }
    }

    smer_stream_o.write(0);
    smer_valid_o.write (0);
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
