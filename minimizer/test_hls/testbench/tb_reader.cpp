
#include <iostream>
#include <hls_stream.h>
#include <ap_int.h>

extern void reader_hls(
    const ap_uint<64>* packed_sequence,
    ap_uint<64> n_bases,
    hls::stream<ap_uint<64>>& base_stream,
    hls::stream<ap_uint<8>>& base_valid
);


int main(){

    const int N = 20;

    uint8_t seq[N] = {
        'A','C','G','T','A','C','G','T',
        'A','C','G','T','A','C','G','T',
        'A','C','G','T'
    };

    ap_uint<64> packed[4];

    for(int i=0;i<4;i++)
        packed[i]=0;

    for(int i=0;i<N;i++){

        int w = i >> 3;
        int b = i & 7;

        packed[w] |= ((ap_uint<64>)seq[i]) << (8*b);
    }

    hls::stream<ap_uint<64>> base_stream;
    hls::stream<ap_uint<8>> base_valid;

    reader_hls(
        packed,
        N,
        base_stream,
        base_valid
    );

    std::cout<<"--- OUTPUT STREAM ---"<<std::endl;

    while(true){

        ap_uint<64> word = base_stream.read();
        ap_uint<8> valid = base_valid.read();

        std::cout<<"valid="<<(int)valid<<std::endl;

        if(valid==0)
            break;

        for(int i=0;i<valid;i++){

            uint8_t c = word.range(8*i+7,8*i);

            std::cout<<(char)c<<" ";
        }

        std::cout<<std::endl;
    }

    return 0;
}
