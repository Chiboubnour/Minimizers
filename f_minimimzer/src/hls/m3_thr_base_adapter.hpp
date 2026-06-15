#pragma once
#include "./header.hpp"
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//  m3 : aligneur de s-mer (template sur PAR_FACTOR), header-only.
//
//  Re-decoupe le flux de bases 2 bits pour aligner les battements sur le 1er
//  s-mer complet : warm-up de (SMER_SIZE-1) bases, livre en first_rounds paquets
//  pleins + 1 paquet de n_last_round bases, puis le reste en paquets pleins.
//
//      first_rounds = (SMER_SIZE-1) / PAR_FACTOR
//      n_last_round = (SMER_SIZE-1) % PAR_FACTOR     (bases du paquet partiel)
//      buff_values  = PAR_FACTOR - n_last_round      (bases reportees, le carry)
//
//  Pour PF=8 : 2, 4, 4 -> cadencage 8,8,4,8,...  ; pour PF=32/64 : first_rounds=0
//  (tout le warm-up tient dans le paquet partiel). Convention : seul "valid" fait
//  foi (donnees des slots invalides = don't-care).
//
template<int PAR_FACTOR>
void m3_thr_base_adapter(
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& base_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid,
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& out_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& out_valid
){
#pragma HLS INLINE off

    constexpr int E            = 2;                              // bits par base
    constexpr int first_rounds = (SMER_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (SMER_SIZE - 1) % PAR_FACTOR;
    constexpr int buff_values  = PAR_FACTOR - n_last_round;
    const ap_uint<PAR_FACTOR> last_valid = ((ap_uint<PAR_FACTOR>)1 << n_last_round) - 1;

    // carry : buff_values bases en attente (donnees + validite)
    ap_uint<E * buff_values> buffer_d = 0;
    ap_uint<    buff_values> buffer_v = 0;

    ap_uint<E * PAR_FACTOR> in_word;
    ap_uint<    PAR_FACTOR> in_valid;

    //
    // first_rounds paquets pleins recopies verbatim
    //
    for (int r = 0; r < first_rounds; r++) {
        in_word  = base_stream.read();
        in_valid = base_valid.read();
        out_stream.write(in_word);
        out_valid.write (in_valid);
    }
    //
    // paquet de warm-up partiel : n_last_round bases valides ; on bufferise les
    // buff_values bases hautes pour l'iteration suivante.
    //
    in_word  = base_stream.read();
    in_valid = base_valid.read();
    out_stream.write(in_word);
    out_valid.write (last_valid);
    buffer_d = in_word.range (E * PAR_FACTOR - 1, E * n_last_round);   // buff_values bases hautes
    buffer_v = in_valid.range(    PAR_FACTOR - 1,     n_last_round);

    //
    // regime permanent : (carry + n_last_round nouvelles) -> paquet plein, et on
    // reporte les buff_values bases hautes du nouveau paquet.
    //
    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        in_word  = base_stream.read();
        in_valid = base_valid.read();

        const ap_uint<E * buff_values>  up_part_d = in_word.range (E * PAR_FACTOR - 1, E * n_last_round);
        const ap_uint<    buff_values>  up_part_v = in_valid.range(    PAR_FACTOR - 1,     n_last_round);
        const ap_uint<E * n_last_round> dw_part_d = in_word.range (E * n_last_round - 1, 0);
        const ap_uint<    n_last_round> dw_part_v = in_valid.range(    n_last_round - 1, 0);

        const ap_uint<E * PAR_FACTOR> to_send_d = (dw_part_d, buffer_d);
        const ap_uint<    PAR_FACTOR> to_send_v = (dw_part_v, buffer_v);

        // terminateur unique : on n'emet que si le paquet porte au moins une base
        if (to_send_v != 0) {
            out_stream.write(to_send_d);
            out_valid.write (to_send_v);
        }

        buffer_d = up_part_d;
        buffer_v = up_part_v;

        if (in_valid == 0) {
            break;
        }
    }

    //
    // paquet terminal
    //
    out_stream.write(0);
    out_valid.write (0);
}
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
