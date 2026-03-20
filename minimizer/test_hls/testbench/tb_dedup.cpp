
#include <iostream>
#include <hls_stream.h>
#include <ap_int.h>

constexpr int SMER = 11;
constexpr int SMER_BITS = 2*SMER;
constexpr int WINDOW = 4;

template<int WINDOW,int SMER_BITS>
void thread_dedup(
    hls::stream<ap_uint<8*SMER_BITS>>& in,
    hls::stream<ap_uint<8>>& valid_i,
    hls::stream<ap_uint<8*SMER_BITS>>& out,
    hls::stream<ap_uint<8>>& valid_o
);

int main() {
    hls::stream<ap_uint<8*SMER_BITS>> in_stream;
    hls::stream<ap_uint<8>> valid_in;

    hls::stream<ap_uint<8*SMER_BITS>> out_stream;
    hls::stream<ap_uint<8>> valid_out;

    // ================================
    // Packet 1
    // ================================
    ap_uint<SMER_BITS> smers1[8] = {
        0x2bb,
        0x2a3,
        0x2a3,
        0x2bb,
        0x2bb,
        0x2a3,
        0x2a3,
        0x2bb
    };

    ap_uint<8*SMER_BITS> packet1 = 0;
    for(int i=0;i<8;i++){
        packet1.range((i+1)*SMER_BITS-1,i*SMER_BITS) = smers1[i];
    }

    in_stream.write(packet1);
    valid_in.write(8);

    // ================================
    // Packet 2
    // ================================
    ap_uint<SMER_BITS> smers2[4] = {
        0x2bb,
        0x2a3,
        0x003,
        0x000
    };

    ap_uint<8*SMER_BITS> packet2 = 0;
    for(int i=0;i<4;i++){
        packet2.range((i+1)*SMER_BITS-1,i*SMER_BITS) = smers2[i];
    }

    in_stream.write(packet2);
    valid_in.write(4);

    // Termination
    in_stream.write(0);
    valid_in.write(0);

    // ================================
    // Call dedup
    // ================================
    thread_dedup<WINDOW,SMER_BITS>(
        in_stream,
        valid_in,
        out_stream,
        valid_out
    );

    // ================================
    // Print results
    // ================================
    std::cout << "\n============================\n";
    std::cout << "DEDUP OUTPUT\n";
    std::cout << "============================\n";

    while(true){
        ap_uint<8*SMER_BITS> packet = out_stream.read();
        ap_uint<8> valid = valid_out.read();

        if(valid == 0) break;

        std::cout << "VALID : " << (int)valid << std::endl;
        std::cout << "PACKET HEX : 0x" << std::hex << packet << std::dec << std::endl;

        for(int i=0;i<valid;i++){
            ap_uint<SMER_BITS> v =
                packet.range((i+1)*SMER_BITS-1, i*SMER_BITS);
            std::cout << "minimizer[" << i << "] = 0x" << std::hex << v << std::dec << std::endl;
        }

        std::cout << std::endl;
    }

    return 0;
}

