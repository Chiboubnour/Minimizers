
#include <iostream>
#include <hls_stream.h>
#include <ap_int.h>

constexpr int S_BASE = 5;
constexpr int S_BITS = 2*S_BASE;

extern void adapter_smer_hls(
    hls::stream<ap_uint<8*S_BITS>>& base_stream,
    hls::stream<ap_uint<8>>&        base_valid,
    hls::stream<ap_uint<8*S_BITS>>& out_stream,
    hls::stream<ap_uint<8>>&        out_valid
);

int main() {
    hls::stream<ap_uint<8*S_BITS>> base_stream;
    hls::stream<ap_uint<8>>        base_valid;
    hls::stream<ap_uint<8*S_BITS>> out_stream;
    hls::stream<ap_uint<8>>        out_valid;


    ap_uint<8*S_BITS> test_data[2];
    ap_uint<8> valid_data[2];

    test_data[0]  = 0x0AEEA3A8EBB; // premier packet
    valid_data[0] = 4;             // 4 smers valides dans ce packet

    test_data[1]  = 0xAEEA3A8EBBAEEA3A8EBB; // deuxième packet
    valid_data[1] = 8;                     // 8 smers valides

    for(int i=0; i<2; i++){
        base_stream.write(test_data[i]);
        base_valid.write(valid_data[i]);
    }

    base_stream.write(0);
    base_valid.write(0);

    adapter_smer_hls(base_stream, base_valid, out_stream, out_valid);

    std::cout << "============================" << std::endl;
    std::cout << "OUTPUT STREAM (adapter_smer_hls)" << std::endl;
    std::cout << "============================" << std::endl;

    while(true){
        ap_uint<8*S_BITS> word = out_stream.read();
        ap_uint<8> valid = out_valid.read();

        if(valid == 0) break;

        std::cout << "VALID : " << (int)valid << std::endl;
        std::cout << "PACKET HEX : 0x" << std::hex << word << std::dec << std::endl;

        for(int i=0; i<valid; i++){
            ap_uint<S_BITS> smer = word.range((i+1)*S_BITS-1, i*S_BITS);
            std::cout << "smer[" << i << "] = 0x" << std::hex << smer << std::dec << std::endl;
        }

        std::cout << std::endl;
    }

    return 0;
}

