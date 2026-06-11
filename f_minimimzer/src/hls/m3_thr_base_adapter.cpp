#include "./m3_thr_base_adapter.hpp"
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
//
//
//
//
void m3_thr_base_adapter(
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& base_stream,    // 2 bits par base
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid,     // 1 bit. par base
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& out_stream,     // 2 bits par base
    hls::stream<ap_uint<     PAR_FACTOR >>& out_valid       // 1 bit. par base
){
#pragma HLS INLINE off
    //constexpr int PAR_FACTOR  =  8;
    //constexpr int SMER_SIZE   = 21;
    constexpr int first_rounds = SMER_SIZE / PAR_FACTOR;                     // 2 rounds       = 16 bases 
    constexpr int n_last_round = SMER_SIZE  - first_rounds * PAR_FACTOR - 1; // 21 - 2 * 8 - 1 =  4 bases
    constexpr int buff_values  = PAR_FACTOR - n_last_round;                  // 8 - 4          =  4 bases

    // buffer 128 bits pour absorber 2 mots de 64 bits (64 bases 2 bits)

    ap_uint<2 * buff_values> buffer_d;
    ap_uint<    buff_values> buffer_v;

    ap_uint<64> in_word ;
    ap_uint<8>  in_valid;

    for(int i = 0; i < first_rounds; i += 1){
        //
        // On lit les entrées du module
        //
        in_word  = base_stream.read();
        in_valid = base_valid.read();
        //
        // On recopie directement les entrée en sortie
        //
        out_stream.write(in_word );
        out_valid.write (in_valid);
        //
    }
    //
    // On lit les données d'entrée
    //
    in_word  = base_stream.read();
    in_valid = base_valid.read();
    //
    // On transmet le dernier paquet vers le module suivante (une partie des données seulement)
    //
    out_stream.write(in_word);
    out_valid.write (0x0F);         // TODO rendre cela générique
    //
    // On bufferise les données pour l'itération suiavnte
    //
    buffer_d = in_word.range (2 * PAR_FACTOR - 1, 2 * buff_values);
    buffer_v = in_valid.range(2 * PAR_FACTOR - 1, 2 * buff_values);;
    //
    //    
    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        //
        in_word  = base_stream.read();
        in_valid = base_valid.read();
        //
        const ap_uint<2 * buff_values> up_part_d = in_word.range (2 * PAR_FACTOR - 1, 2 * buff_values);
        const ap_uint<    buff_values> up_part_v = in_valid.range(2 * PAR_FACTOR - 1, 2 * buff_values);
        //
        const ap_uint<2 * buff_values> dw_part_d = in_word.range (2 * buff_values - 1, 0);
        const ap_uint<    buff_values> dw_part_v = in_valid.range(2 * buff_values - 1, 0);
        //
        const ap_uint<2 * buff_values> to_send_d = (dw_part_d, buffer_d);
        const ap_uint<    buff_values> to_send_v = (dw_part_v, buffer_v);
        //
        out_stream.write(to_send_d);
        out_valid.write (to_send_v);
        //
        buffer_d = up_part_d;   // on memorise les données restante pour l'itération suivante
        buffer_v = up_part_v;   // les data et les bits de validation
        //
        if (in_valid == 0)
        {
            break;
        }
    }
    //
    out_stream.write(0);
    out_valid.write (0);
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