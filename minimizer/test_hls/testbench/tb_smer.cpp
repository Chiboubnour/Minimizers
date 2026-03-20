#include <iostream>
#include <hls_stream.h>
#include <ap_int.h>

#define S_BASE 5

extern void smer_thread_hls(
    hls::stream<ap_uint<64>>&,
    hls::stream<ap_uint<8>>&,
    hls::stream<ap_uint<16*S_BASE>>&,
    hls::stream<ap_uint<8>>&
);

extern void smer_generate_thread(
    hls::stream<ap_uint<64>>&,
    hls::stream<ap_uint<8>>&,
    hls::stream<ap_uint<16*S_BASE>>&,
    hls::stream<ap_uint<8>>&
);

extern void hash_thread(
    hls::stream<ap_uint<16*S_BASE>>&,
    hls::stream<ap_uint<8>>&,
    hls::stream<ap_uint<16*S_BASE>>&,
    hls::stream<ap_uint<8>>&
);

int main()
{
    const char seq[] = "ACGTACGTACGTACGT";
    const int N = sizeof(seq)-1;

    std::cout<<"Sequence: "<<seq<<"\n\n";

    // duplication streams
    hls::stream<ap_uint<64>> base_stream1, base_stream2;
    hls::stream<ap_uint<8>>  base_valid1, base_valid2;

    ap_uint<64> word = 0;

    for(int i=0;i<N;i++){
        word |= ((ap_uint<64>)seq[i]) << (8*(i%8));

        if(i%8==7){
            base_stream1.write(word);
            base_valid1.write(8);

            base_stream2.write(word);
            base_valid2.write(8);

            word = 0;
        }
    }

    int rem = N % 8;

    if(rem!=0){
        base_stream1.write(word);
        base_valid1.write(rem);

        base_stream2.write(word);
        base_valid2.write(rem);
    }

    // end
    base_stream1.write(0); base_valid1.write(0);
    base_stream2.write(0); base_valid2.write(0);

    // outputs
    hls::stream<ap_uint<16*S_BASE>> out1;
    hls::stream<ap_uint<8>> valid1;

    hls::stream<ap_uint<16*S_BASE>> smer_stream;
    hls::stream<ap_uint<8>> smer_valid;

    hls::stream<ap_uint<16*S_BASE>> out2;
    hls::stream<ap_uint<8>> valid2;

    // run version 1
    smer_thread_hls(base_stream1, base_valid1, out1, valid1);

    // run version 2
    smer_generate_thread(base_stream2, base_valid2, smer_stream, smer_valid);
    hash_thread(smer_stream, smer_valid, out2, valid2);

    // compare
    std::cout<<"========================\n";
    std::cout<<"COMPARAISON\n";
    std::cout<<"========================\n";

    int pkt = 0;

    while(true){
        auto p1 = out1.read();
        auto v1 = valid1.read();

        auto p2 = out2.read();
        auto v2 = valid2.read();

        std::cout<<"Packet "<<pkt++<<"\n";
        std::cout<<"V1="<<(int)v1<<" V2="<<(int)v2<<"\n";

        if(v1 != v2){
            std::cout<<"❌ VALID mismatch\n";
        }

        for(int i=0;i<8;i++){
            if(i < v1 && i < v2){
                auto s1 = p1.range((i+1)*2*S_BASE-1,i*2*S_BASE);
                auto s2 = p2.range((i+1)*2*S_BASE-1,i*2*S_BASE);

                if(s1 != s2){
                    std::cout<<"❌ mismatch @ "<<i
                             <<" : "<<std::hex<<s1
                             <<" vs "<<s2<<std::dec<<"\n";
                }
            }
        }

        if(v1==0 && v2==0)
            break;
    }

    std::cout<<"FIN\n";
    return 0;
}
