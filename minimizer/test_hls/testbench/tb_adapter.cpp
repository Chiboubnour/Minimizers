
#include <iostream>
#include <hls_stream.h>
#include <ap_int.h>

extern void adapter_hls(
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>&  base_valid,
    hls::stream<ap_uint<64>>& out_stream,
    hls::stream<ap_uint<8>>&  out_valid
);

int main(){

    hls::stream<ap_uint<64>> base_stream;
    hls::stream<ap_uint<8>>  base_valid;

    hls::stream<ap_uint<64>> out_stream;
    hls::stream<ap_uint<8>>  out_valid;

    // =============================
    // Données d'entrée (simulation)
    // =============================

    // word 1 : 8 bases
    ap_uint<64> w1 = 0;
    w1.range(7,0)   = 'A';
    w1.range(15,8)  = 'C';
    w1.range(23,16) = 'G';
    w1.range(31,24) = 'T';
    w1.range(39,32) = 'A';
    w1.range(47,40) = 'C';
    w1.range(55,48) = 'G';
    w1.range(63,56) = 'T';

    // word 2 : 8 bases
    ap_uint<64> w2 = w1;

    // word 3 : 4 bases
    ap_uint<64> w3 = 0;
    w3.range(7,0)   = 'A';
    w3.range(15,8)  = 'C';
    w3.range(23,16) = 'G';
    w3.range(31,24) = 'T';

    // écriture dans le stream
    base_stream.write(w1);
    base_valid.write(8);

    base_stream.write(w2);
    base_valid.write(8);

    base_stream.write(w3);
    base_valid.write(4);

    // fin du stream
    base_stream.write(0);
    base_valid.write(0);


    adapter_hls(
        base_stream,
        base_valid,
        out_stream,
        out_valid
    );

    std::cout << "--- OUTPUT ---" << std::endl;

    while(true){

        ap_uint<64> word = out_stream.read();
        ap_uint<8> valid = out_valid.read();

        std::cout << "valid=" << (int)valid << std::endl;

        if(valid == 0)
            break;

        for(int i=0;i<valid;i++){

            uint8_t c = word.range(8*i+7,8*i);

            std::cout << (char)c << " ";
        }

        std::cout << std::endl;
    }

    return 0;
}

