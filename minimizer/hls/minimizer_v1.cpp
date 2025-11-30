#include <ap_int.h>
 #include <hls_stream.h>
 #include <cstdint>
 
 ///////////////////////////////////////////////////////////////////////////////////
 // Paramètres
 ///////////////////////////////////////////////////////////////////////////////////
 #define SMER_SIZE 56
 #define WINDOW_SIZE 16
 #define DATA_DEPTH 1024
 #define MEM_UNIT 64
 
inline ap_uint<2> nucl_encode_v1(const char nucl)
 {
     return (nucl >> 1) & 0b11;
 }

template<int width>
inline ap_uint<64>min_v1(const ap_uint<width> a, const ap_uint<width> b)
{
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE
#endif
    return (a < b) ? a : b;
}

 
 inline ap_uint<64> mask_right(int numbits) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE
#endif
     return (numbits >= MEM_UNIT) ? ~0ULL : ((1ULL << numbits) - 1ULL);
 }
 

template <int resu_size>
inline ap_uint<resu_size> hash_u64(uint64_t key) 
 {
     key = (~key + (key << 21)) /*& mask*/; // key = (key << 21) - key - 1;
     key ^= key >> 24;
     key = ((key + (key << 3)) + (key << 8)) /*& mask*/; // key * 265
     key ^= key >> 14;
     key = ((key + (key << 2)) + (key << 4)) /*& mask*/; // key * 21
     key ^= key >> 28;
     key = (key + (key << 31)) /*& mask*/;
     return key;
 }
 #if 0
 void thread_reader_v2(
     const ap_uint<64>* packed_sequence,
     ap_uint<64> n_bases,
     hls::stream< ap_uint<24> >& stream_o
 ) {
     const int n_words = (int)((n_bases + 7) / 8);
 
     for (int i = 0; i < n_words; ++i) {
#ifdef _PRAGMA_HLS_
     #pragma HLS PIPELINE II=1
#endif
         const ap_uint<64> word_8b = packed_sequence[i];
 
         ap_uint<24> word_3b = 0;
         bool all_valid = true;
 
         // 8 caractères ASCII -> 8 triplets (2 bits encodés + 1 bit valid)
         for (int j = 0; j < 8; ++j) {
#ifdef _PRAGMA_HLS_
         #pragma HLS UNROLL
#endif
             const int idx = i * 8 + j;
             ap_uint<1> valid = (idx < (int)n_bases) ? ap_uint<1>(1) : ap_uint<1>(0);
             ap_uint<8> c     = valid ? word_8b.range(8*(j+1)-1, 8*j) : ap_uint<8>(0);
 
             const ap_uint<2> enc = (c >> 1) & 0x3;
             word_3b.range(3*j+1, 3*j) = enc;  // 2 bits
             word_3b[3*j+2]            = valid;
 
             all_valid &= (bool)valid;
         }
 
         stream_o.write(word_3b);
          if (!all_valid) break;
     }
 
     if ((n_bases & 7) == 0) {
#ifdef _PRAGMA_HLS_
        #pragma HLS PIPELINE II=1
#endif
            stream_o.write( (ap_uint<24>)0 );
     }
 }
 #else
 template <int smer_size = 19>
void thread_reader(
    const uint64_t* packed_sequence,
    hls::stream<uint8_t>& stream_o
) {
#ifdef _PRAGMA_HLS_
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1
#endif
    bool stop = false;
    for (int i = 0; i < 64 && !stop; i++) {
        uint64_t word = packed_sequence[i];
        for (int j = 0; j < 8; j++) {
            char c = word & 0xFF;
            if (c == 0) {
                stop = true;
                break;
            }
            stream_o.write(c);
            word >>= 8;
        }
    }
    // Terminate the stream
    stream_o.write(0x00);
}
 #endif
 void thread_smer_v2(
     hls::stream< ap_uint<24> >& stream_i,
     ap_uint<64> n_bases,
     hls::stream< ap_uint<SMER_SIZE> >& stream_o
 ) {
     constexpr int smer = SMER_SIZE / 2;
     const ap_uint<64> HASH_MASK = mask_right(SMER_SIZE);
 
     ap_uint<SMER_SIZE> current_smer = 0;
     ap_uint<SMER_SIZE> cur_inv_smer = 0;
     ap_uint<24> word_24b = 0;
 
     for (int i = 0; i < smer - 1; i++) {
#ifdef _PRAGMA_HLS_
        #pragma HLS PIPELINE II=1
#endif       
    if ((i & 0x07) == 0){
             word_24b = stream_i.read();
         }
         const ap_uint<2> c_nucl = word_24b.range(1,0);
 
         current_smer <<= 2;
         current_smer(1,0) = c_nucl;
 
         cur_inv_smer >>= 2;
         cur_inv_smer(SMER_SIZE-1, SMER_SIZE-2) = (0x2 ^ c_nucl);
 
         word_24b >>= 3;
     }
 
     const int last_pos = (int)n_bases - 1;
     for (int i = smer - 1; i <= last_pos; i++) {
        
#ifdef _PRAGMA_HLS_
        #pragma HLS PIPELINE II=1
#endif
        if ((i & 0x07) == 0) word_24b = stream_i.read();
 
         const ap_uint<2> c_nucl = word_24b.range(1,0);
         const bool valid = word_24b[2];
 
         current_smer <<= 2;
         current_smer(1,0) = c_nucl;
 
         cur_inv_smer = (cur_inv_smer >> 2) | ((ap_uint<SMER_SIZE>)((0x2 ^ c_nucl)) << (SMER_SIZE-2));
 
         const ap_uint<SMER_SIZE> vmin  = min_v1(current_smer, cur_inv_smer);
         const ap_uint<SMER_SIZE> vhash = hash_u64(vmin, HASH_MASK);
 
         word_24b >>= 3;
 
         if (valid) {
             stream_o.write(vhash);
         } else {
             stream_o.write(0x00);
             return;
         }
     }
 
     stream_o.write(0x00);
 }
 
 
 ///////////////////////////////////////////////////////////////////////////////////
 // Dedup within window (tree reduction pour minimum)
 ///////////////////////////////////////////////////////////////////////////////////
 void thread_dedup_v2(
     hls::stream< ap_uint<SMER_SIZE> >& stream_i,
     hls::stream< ap_uint<SMER_SIZE> >& stream_o
 ) {
     ap_uint<SMER_SIZE> buffer[WINDOW_SIZE];
     for (int p = 0; p < WINDOW_SIZE; p++) {
     #pragma HLS PIPELINE II=1
         buffer[p] = stream_i.read();
     }
 
     ap_uint<SMER_SIZE> lastElement = (ap_uint<SMER_SIZE>)(-1);
 
     while (true) {
     #pragma HLS PIPELINE II=1
         const ap_uint<SMER_SIZE> vhash = stream_i.read();
         if (vhash == 0) { stream_o.write(0x00); break; }
 
         ap_uint<SMER_SIZE> minz = vhash;
         for (int p = 0; p < WINDOW_SIZE; p++) {  
         #pragma HLS UNROLL
             minz = (buffer[p] < minz) ? buffer[p] : minz;
         }
         for (int p = 0; p < WINDOW_SIZE - 1; p++) {
         #pragma HLS UNROLL
             buffer[p] = buffer[p+1];
         }
         buffer[WINDOW_SIZE-1] = vhash;
 
         if (lastElement != minz) {
             stream_o.write(minz);
             lastElement = minz;
         }
     }
 }
 
 ///////////////////////////////////////////////////////////////////////////////////
 // Store
 ///////////////////////////////////////////////////////////////////////////////////
 void thread_store_v2(
     hls::stream< ap_uint<SMER_SIZE> >& stream_i,
     ap_uint<64>* tab_hash,
     ap_uint<64>* nElements
 ) {
     int cnt = 0;
     while (true) {
     #pragma HLS PIPELINE II=1
         const ap_uint<SMER_SIZE> vHash = stream_i.read();
         if (vHash == 0) break;
         tab_hash[cnt++] = vHash;
     }
     *nElements = cnt;
 }
 
 ///////////////////////////////////////////////////////////////////////////////////
 // Top-level
 ///////////////////////////////////////////////////////////////////////////////////
 #define SMER 19
 
 void minimizer(
     const ap_uint<64>* packed_sequence,
     ap_uint<64> n,
     ap_uint<64>* tab_hash,
     ap_uint<64>* nMinizrs
 ) {
#ifdef _PRAGMA_HLS_
     #pragma HLS INTERFACE mode=m_axi port=packed_sequence
     #pragma HLS INTERFACE mode=m_axi port=tab_hash
     #pragma HLS INTERFACE mode=s_axilite port=nMinizrs
     #pragma HLS INTERFACE  mode=s_axilite port=n
     #pragma HLS INTERFACE mode=s_axilite port=return bundle=control
     #pragma HLS DATAFLOW
#endif
 
 
     hls::stream< ap_uint<       24>, DATA_DEPTH > fifo_1;
     hls::stream< ap_uint<SMER_SIZE>, DATA_DEPTH > fifo_2;
     hls::stream< ap_uint<SMER_SIZE>, DATA_DEPTH > fifo_3;
 
     thread_reader_v2(packed_sequence, n, fifo_1);
     thread_smer_v2(fifo_1, n, fifo_2);
     thread_dedup_v2(fifo_2, fifo_3);
     thread_store_v2(fifo_3, tab_hash, nMinizrs);
 }