#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include "minimizer.hpp"


template<int width>
inline ap_uint<64> min_v1(const ap_uint<width> a, const ap_uint<width> b)
{
#pragma HLS INLINE
    return (a < b) ? a : b;
}

template <int resu_size>
ap_uint<resu_size> hash_u64(ap_uint<64> key)
{
#pragma HLS INLINE off
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return key;
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
template<int PAR_FACTOR>
inline ap_uint<PAR_FACTOR> valid_bits(const ap_uint<64> reste)
{
    if (reste >= PAR_FACTOR) return ~ap_uint<PAR_FACTOR>(0);
    return (ap_uint<PAR_FACTOR>(1) << reste) - 1;
}


template<int PAR_FACTOR>
void thr_reader_mem(
    const ap_uint<MEM_WIDTH>*          base_ptr_i,
    const ap_uint<64>                  n_bases_i,
    hls::stream<ap_uint<MEM_WIDTH>>&   mem_stream_o
){
#pragma HLS INLINE off
    const ap_uint<64> n_words = (n_bases_i + PAR_FACTOR - 1) / PAR_FACTOR;
    const ap_uint<64> n_beats = (n_words + MEM_RATIO - 1) / MEM_RATIO;

    for (ap_uint<64> b = 0; b < n_beats; b++)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        mem_stream_o.write(base_ptr_i[b]);
    }
}

// thr_unpack : decoupe chaque beat de MEM_WIDTH bits en MEM_RATIO mots de
// 8*PAR_FACTOR bits
template<int PAR_FACTOR>
void thr_unpack(
    hls::stream<ap_uint<MEM_WIDTH>>&          mem_stream_i,
    const ap_uint<64>                         n_bases_i,
    hls::stream<ap_uint<8 * PAR_FACTOR>>&     base_stream_o,
    hls::stream<ap_uint<PAR_FACTOR>>&         base_valid_o
){
#pragma HLS INLINE off
    const ap_uint<64> n_words = (n_bases_i + PAR_FACTOR - 1) / PAR_FACTOR;
    ap_uint<64> base_cnt = n_bases_i;
    ap_uint<MEM_WIDTH> beat = 0;

    for (ap_uint<64> w = 0; w < n_words; w++)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<64> sub = w % MEM_RATIO;
        if (sub == 0) {
            beat = mem_stream_i.read();
        }

        const ap_uint<8 * PAR_FACTOR> word  = beat.range((sub + 1) * 8 * PAR_FACTOR - 1, sub * 8 * PAR_FACTOR);
        const ap_uint<     PAR_FACTOR> valid = valid_bits<PAR_FACTOR>(base_cnt);

        base_stream_o.write(word);
        base_valid_o.write(valid);

        base_cnt -= PAR_FACTOR;
    }
    base_stream_o.write(0);
    base_valid_o.write(0);
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
template<int PAR_FACTOR>
void thr_adapter_hls(
    hls::stream<ap_uint< 8 * PAR_FACTOR >>& base_stream_i,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid_i,
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& base_stream_o,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid_o
){
#pragma HLS INLINE off
    ap_uint<2 * PAR_FACTOR> buffer = 0;
    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        const ap_uint<8 * PAR_FACTOR> in_word  = base_stream_i.read();
        const ap_uint<    PAR_FACTOR> in_valid = base_valid_i.read();
        for (int i = 0; i < PAR_FACTOR; i++)
        {
            const ap_uint<8> nucl  = in_word.range(8 * i + 7, 8 * i);
            const ap_uint<2> enc2b = (nucl >> 1) & 3;
            buffer.range(2 * i + 1, 2 * i) = enc2b;
        }
        base_stream_o.write(buffer);
        base_valid_o.write(in_valid);
        if (in_valid == 0) break;
    }
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
template<int PAR_FACTOR>
void m3_thr_base_adapter(
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& base_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid,
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& out_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& out_valid
){
#pragma HLS INLINE off
    constexpr int E            = 2;
    constexpr int SHIFT        = SMER_SIZE - 1;
    constexpr int first_rounds = SHIFT / PAR_FACTOR;
    constexpr int offset       = SHIFT % PAR_FACTOR;

    ap_uint<E * PAR_FACTOR> buffer_d = 0;
    ap_uint<    PAR_FACTOR> buffer_v = 0;

    for (int r = 0; r < first_rounds; r++) {
        auto d = base_stream.read();
        auto v = base_valid.read();
        out_stream.write(d);
        out_valid.write(v);
    }

    {
        auto d = base_stream.read();
        auto v = base_valid.read();

        if constexpr (offset == 0) {
            out_stream.write(d);
            out_valid.write(v);
        } else {
            ap_uint<PAR_FACTOR> first_valid = 0;
            for (int i = 0; i < offset; i++)
                first_valid[i] = 1;

            out_stream.write(d);
            out_valid.write(first_valid);

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
#pragma HLS UNROLL
                buffer_d.range(E*(i+1)-1, E*i) = d.range(E*(offset+i+1)-1, E*(offset+i));
                buffer_v[i] = v[offset + i];
            }
        }
    }

    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        auto d = base_stream.read();
        auto v = base_valid.read();

        if (v == 0) break;

        if constexpr (offset == 0) {
            out_stream.write(d);
            out_valid.write(v);
        } else {
            ap_uint<E * PAR_FACTOR> out_data = 0;
            ap_uint<    PAR_FACTOR> out_v    = 0;

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                out_data.range(E*(i+1)-1, E*i) = buffer_d.range(E*(i+1)-1, E*i);
                out_v[i] = buffer_v[i];
            }
            for (int i = 0; i < offset; i++) {
                out_data.range(E*(PAR_FACTOR-offset+i+1)-1, E*(PAR_FACTOR-offset+i)) = d.range(E*(i+1)-1, E*i);
                out_v[PAR_FACTOR-offset+i] = v[i];
            }

            if (out_v != 0) {
                out_stream.write(out_data);
                out_valid.write(out_v);
            }

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                buffer_d.range(E*(i+1)-1, E*i) = d.range(E*(offset+i+1)-1, E*(offset+i));
                buffer_v[i] = v[offset + i];
            }
        }
    }

    out_stream.write(0);
    out_valid.write(0);
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
template<int PAR_FACTOR>
void thr_smer_gen(
    hls::stream<ap_uint<2 *             PAR_FACTOR>>& base_stream_i,
    hls::stream<ap_uint<                PAR_FACTOR>>& base_valid_i,
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_valid_o
){
#pragma HLS INLINE off
    constexpr int SMER_BITS    = 2 * SMER_SIZE;
    constexpr int HIST_BITS    = 2 * (SMER_SIZE - 1);
    constexpr int WIN_BITS     = HIST_BITS + 2 * PAR_FACTOR;
    constexpr int first_rounds = (SMER_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (SMER_SIZE - 1) % PAR_FACTOR;
    ap_uint<HIST_BITS> memory = 0;
    ap_uint<2 * PAR_FACTOR> in_word;
    ap_uint<    PAR_FACTOR> in_valid;
    for (int r = 0; r < first_rounds; r++) {
        in_word  = base_stream_i.read();
        in_valid = base_valid_i.read();
        for (int i = 0; i < PAR_FACTOR; i++) {
            memory = (memory << 2);
            memory(1, 0) = in_word.range(2 * i + 1, 2 * i);
        }
    }
    if constexpr (n_last_round != 0) {
        in_word  = base_stream_i.read();
        in_valid = base_valid_i.read();
        for (int i = 0; i < n_last_round; i++) {
            memory = (memory << 2);
            memory(1, 0) = in_word.range(2 * i + 1, 2 * i);
        }
    }
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        in_word  = base_stream_i.read();
        in_valid = base_valid_i.read();
        if (in_valid == 0) break;
        ap_uint<2 * PAR_FACTOR> new_rev = 0;
        for (int j = 0; j < PAR_FACTOR; j++) {
            new_rev.range(2*(PAR_FACTOR-1-j)+1, 2*(PAR_FACTOR-1-j)) =
                in_word.range(2*j+1, 2*j);
        }
        const ap_uint<WIN_BITS> win = (memory, new_rev);
        ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_out = 0;
        for (int i = 0; i < PAR_FACTOR; i++)
        {
            const ap_uint<SMER_BITS> fwd = win.range( 2*(PAR_FACTOR-1-i)+2*SMER_SIZE-1, 2*(PAR_FACTOR-1-i));
            ap_uint<SMER_BITS> rev = 0;
            for (int t = 0; t < SMER_SIZE; t++)
            {
                const ap_uint<2> base = fwd.range(2*t+1, 2*t);
                rev.range(2*(SMER_SIZE-1-t)+1, 2*(SMER_SIZE-1-t)) =  base ^ ap_uint<2>(0x2);
            }
            const ap_uint<SMER_BITS> canon = (fwd < rev) ? fwd : rev;
            packed_out.range((i+1)*SMER_BITS-1, i*SMER_BITS) = canon;
        }
        smer_stream_o.write(packed_out);
        smer_valid_o.write(in_valid);
        memory = win.range(HIST_BITS - 1, 0);
    }
    smer_stream_o.write(0);
    smer_valid_o.write(0);
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
template<int PAR_FACTOR>
void thr_hash(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  hash_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_valid_o
){
#pragma HLS INLINE off
    constexpr int S_BITS = 2 * SMER_SIZE;
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        const ap_uint<2*SMER_SIZE*PAR_FACTOR> packed_smer = smer_stream_i.read();
        const ap_uint<            PAR_FACTOR> valid       = smer_valid_i.read();
        if (valid == 0) {
            hash_stream_o.write(0);
            hash_valid_o.write(0);
            break;
        }
        ap_uint<2*SMER_SIZE*PAR_FACTOR> packed_hash = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            const ap_uint<S_BITS> smer = packed_smer.range((i+1)*S_BITS-1, i*S_BITS);
            const ap_uint<S_BITS> h    = hash_u64<S_BITS>((ap_uint<64>)smer);
            packed_hash.range((i+1)*S_BITS-1, i*S_BITS) = h;
        }
        hash_stream_o.write(packed_hash);
        hash_valid_o.write(valid);
    }
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
template<int PAR_FACTOR>
void thr_adapter_smer(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_o
){
#pragma HLS INLINE off
    constexpr int S_BITS       = 2 * SMER_SIZE;
    constexpr int SHIFT        = WINDOW_SIZE - 1;
    constexpr int first_rounds = SHIFT / PAR_FACTOR;
    constexpr int offset       = SHIFT % PAR_FACTOR;

    ap_uint<S_BITS * PAR_FACTOR> buffer_d = 0;
    ap_uint<PAR_FACTOR>          buffer_v = 0;

    for (int r = 0; r < first_rounds; r++) {
        auto d = smer_stream_i.read();
        auto v = smer_valid_i.read();
        smer_stream_o.write(d);
        smer_valid_o.write(v);
    }

    {
        auto d = smer_stream_i.read();
        auto v = smer_valid_i.read();

        if constexpr (offset == 0) {
            smer_stream_o.write(d);
            smer_valid_o.write(v);
        } else {
            ap_uint<PAR_FACTOR> first_valid = 0;
            for (int i = 0; i < offset; i++)
                first_valid[i] = 1;

            smer_stream_o.write(d);
            smer_valid_o.write(first_valid);

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                buffer_d.range((i+1)*S_BITS-1, i*S_BITS) =
                    d.range((offset+i+1)*S_BITS-1, (offset+i)*S_BITS);
                buffer_v[i] = v[offset + i];
            }
        }
    }

    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        auto d = smer_stream_i.read();
        auto v = smer_valid_i.read();

        if (v == 0) break;

        if constexpr (offset == 0) {
            smer_stream_o.write(d);
            smer_valid_o.write(v);
        } else {
            ap_uint<S_BITS*PAR_FACTOR> out_data = 0;
            ap_uint<PAR_FACTOR>        out_v    = 0;

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                out_data.range((i+1)*S_BITS-1, i*S_BITS) = buffer_d.range((i+1)*S_BITS-1, i*S_BITS);
                out_v[i] = buffer_v[i];
            }
            for (int i = 0; i < offset; i++) {
                out_data.range((PAR_FACTOR-offset+i+1)*S_BITS-1, (PAR_FACTOR-offset+i)*S_BITS) =   d.range((i+1)*S_BITS-1, i*S_BITS);
                out_v[PAR_FACTOR-offset+i] = v[i];
            }

            if (out_v != 0) {
                smer_stream_o.write(out_data);
                smer_valid_o.write(out_v);
            }

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                buffer_d.range((i+1)*S_BITS-1, i*S_BITS) =   d.range((offset+i+1)*S_BITS-1, (offset+i)*S_BITS);
                buffer_v[i] = v[offset + i];
            }
        }
    }

    smer_stream_o.write(0);
    smer_valid_o.write(0);
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
template<int PAR_FACTOR>
void thr_min_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  hash_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_o
){
#pragma HLS INLINE off
    constexpr int S_BITS       = 2 * SMER_SIZE;
    constexpr int HIST_BITS    = (WINDOW_SIZE - 1) * S_BITS;
    constexpr int WIN_BITS     = HIST_BITS + PAR_FACTOR * S_BITS;
    constexpr int first_rounds = (WINDOW_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (WINDOW_SIZE - 1) % PAR_FACTOR;
    ap_uint<HIST_BITS> mem = 0;
    ap_uint<S_BITS * PAR_FACTOR> in_word;
    ap_uint<         PAR_FACTOR> in_valid;
    for (int r = 0; r < first_rounds; r++) {
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        mem.range(S_BITS*PAR_FACTOR*(r+1)-1, S_BITS*PAR_FACTOR*r) = in_word;
    }
    if constexpr (n_last_round != 0) {
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        mem.range(HIST_BITS-1, S_BITS*PAR_FACTOR*first_rounds) =
            in_word.range(n_last_round*S_BITS-1, 0);
    }
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        if (in_valid == 0) break;
        const ap_uint<WIN_BITS> win = (in_word, mem);
        ap_uint<S_BITS * PAR_FACTOR> packed_out = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            ap_uint<S_BITS> m = win.range((i+1)*S_BITS-1, i*S_BITS);
            for (int w = 1; w < WINDOW_SIZE; w++) {
                const ap_uint<S_BITS> v =  win.range((i+w+1)*S_BITS-1, (i+w)*S_BITS);
                if (v < m) m = v;
            }
            packed_out.range((i+1)*S_BITS-1, i*S_BITS) = m;
        }
        min_stream_o.write(packed_out);
        min_valid_o.write(in_valid);
        mem = win.range(WIN_BITS-1, PAR_FACTOR*S_BITS);
    }
    min_stream_o.write(0);
    min_valid_o.write(0);
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

template<int PAR_FACTOR>
void thr_dedup_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  dedup_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  dedup_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * SMER_SIZE;   // un minimizer

    ap_uint<S_BITS> last     = 0;
    bool            has_last = false;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<2 * SMER_SIZE * PAR_FACTOR> in_word  = min_stream_i.read();
        const ap_uint<                PAR_FACTOR> in_valid = min_valid_i.read();

        if (in_valid == 0) break;


        ap_uint<S_BITS> scan_val[PAR_FACTOR];
        bool            scan_has[PAR_FACTOR];

        for (int i = 0; i < PAR_FACTOR; i++) {
            scan_val[i] = in_word.range((i + 1) * S_BITS - 1, i * S_BITS);
            scan_has[i] = in_valid[i];
        }

        for (int d = 1; d < PAR_FACTOR; d <<= 1) {
            ap_uint<S_BITS> nxt_val[PAR_FACTOR];
            bool            nxt_has[PAR_FACTOR];

            for (int i = 0; i < PAR_FACTOR; i++) {
                const int  j    = (i >= d) ? (i - d) : i;      // toujours dans [0,PAR_FACTOR)
                const bool take = (i >= d) && !scan_has[i];
                nxt_val[i] = take ? scan_val[j] : scan_val[i];
                nxt_has[i] = take ? scan_has[j] : scan_has[i];
            }
            for (int i = 0; i < PAR_FACTOR; i++) {
                scan_val[i] = nxt_val[i];
                scan_has[i] = nxt_has[i];
            }
        }

        ap_uint<PAR_FACTOR> is_new = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (in_valid[i]) {
                const bool            prev_has = (i == 0) ? has_last : (scan_has[i - 1] || has_last);
                const ap_uint<S_BITS> prev_val = (i == 0) ? last     : (scan_has[i - 1] ? scan_val[i - 1] : last);
                const ap_uint<S_BITS> v        = in_word.range((i + 1) * S_BITS - 1, i * S_BITS);
                if (!prev_has || v != prev_val)
                    is_new[i] = 1;
            }
        }

        // Densite theorique en sortie ~ 2/(WINDOW_SIZE+1) : la plupart des
        // cycles n'ont aucune voie retenue, d'ou l'ecriture conditionnelle
        // (evite de saturer le stream aval avec des mots tout a zero).
        if (is_new != 0) {
            dedup_stream_o.write(in_word);
            dedup_valid_o.write(is_new);
        }

        if (scan_has[PAR_FACTOR - 1]) {
            last     = scan_val[PAR_FACTOR - 1];
            has_last = true;
        }
    }

    dedup_stream_o.write(0);
    dedup_valid_o.write(0);
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
template<int PAR_FACTOR>
void thr_compact(
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& in_stream,
    hls::stream<ap_uint<PAR_FACTOR>>&                 in_valid,
    hls::stream<ap_uint<64 * PAR_FACTOR>>&            out_stream,
    hls::stream<ap_uint<8>>&                          out_count)
{
#pragma HLS INLINE off
    constexpr int S_BITS = 2 * SMER_SIZE;

COMPACT_LOOP:
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        ap_uint<S_BITS * PAR_FACTOR> word  = in_stream.read();
        ap_uint<PAR_FACTOR>          valid = in_valid.read();

        if (valid == 0) {
            out_stream.write(0);
            out_count.write(0);
            break;
        }


        ap_uint<4> scan[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++) {
            scan[i] = valid[i];
        }
        for (int d = 1; d < PAR_FACTOR; d <<= 1) {
            ap_uint<4> nxt[PAR_FACTOR];
            for (int i = 0; i < PAR_FACTOR; i++) {
                nxt[i] = (i >= d) ? (ap_uint<4>)(scan[i] + scan[i - d]) : scan[i];
            }
            for (int i = 0; i < PAR_FACTOR; i++) {
                scan[i] = nxt[i];
            }
        }

        const ap_uint<4> count = scan[PAR_FACTOR - 1];


        ap_uint<64 * PAR_FACTOR> out_word = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            const ap_uint<4>  pos  = scan[i] - valid[i];   // prefix exclusif
            const ap_uint<64> lane = word.range((i + 1) * S_BITS - 1, i * S_BITS);
            if (valid[i]) {
                out_word |= (ap_uint<64 * PAR_FACTOR>(lane) << (64 * pos));
            }
        }

        out_stream.write(out_word);
        out_count.write((ap_uint<8>)count);
    }
}
//
//
//
//

//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
// thr_flatten : convertit le format compacte de thr_compact (paquet
// jusqu'a PAR_FACTOR elements/cycle + compteur) en un flux SCALAIRE
// valeur+enable, 1 element/cycle, terminé par un enable=false --
//
template<int PAR_FACTOR>
void thr_flatten(
    hls::stream<ap_uint<64 * PAR_FACTOR>>& elem_stream_i,
    hls::stream<ap_uint<8>>&               elem_count_i,
    hls::stream<ap_uint<GROUP_BITS>>&      out_v,
    hls::stream<bool>&                     out_e,
    hls::stream<ap_uint<64>>&              total_o)
{
#pragma HLS INLINE off

    ap_uint<64 * PAR_FACTOR> word   = 0;
    ap_uint<8>                remain = 0;
    ap_uint<8>                pos    = 0;

    ap_uint<GROUP_BITS> group_buf = 0;
    ap_uint<8>          group_pos = 0;   // 0..GROUP_W-1
    ap_uint<64>         total     = 0;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        if (remain == 0) {
            word = elem_stream_i.read();
            const ap_uint<8> count = elem_count_i.read();
            if (count == 0) {
                if (group_pos != 0) {

                    out_v.write(group_buf);
                    out_e.write(true);
                }
                out_v.write(0);
                out_e.write(false);   // sentinelle de fin
                total_o.write(total);
                break;
            }
            remain = count;
            pos    = 0;
        }

        const ap_uint<64> elem = word.range(64 * (pos + 1) - 1, 64 * pos);
        pos++;
        remain--;
        total++;

        group_buf.range(64 * (group_pos + 1) - 1, 64 * group_pos) = elem;
        group_pos++;

        if (group_pos == GROUP_W) {
            out_v.write(group_buf);
            out_e.write(true);
            group_buf = 0;
            group_pos = 0;
        }
    }
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
constexpr int BURST_LEN = 512;
void thr_write_burst(
    hls::stream<ap_uint<GROUP_BITS>>& inv,
    hls::stream<bool>&                ine,
    ap_uint<GROUP_BITS>*               dst,
    ap_uint<64>                        n_bases)
{
#pragma HLS INLINE off

    bool                lastE = false;
    ap_uint<GROUP_BITS> x     = 0;
    bool                e     = false;

    const ap_uint<64> n_groups_max = (n_bases + GROUP_W - 1) / GROUP_W;

    for (ap_uint<64> i = 0; i < n_groups_max; i += BURST_LEN)
    {
#pragma HLS LOOP_TRIPCOUNT min=64 max=64 avg=64
#pragma HLS DEPENDENCE variable=dst class=array inter false
        for (int j = 0; j < BURST_LEN; j++)
        {
#pragma HLS PIPELINE II=1
            if (!lastE) {
                x = inv.read();
                e = ine.read();
                if (!e) lastE = true;
            }
            dst[i + j] = x;   // ecriture inconditionnelle (hors du if) --
                               // necessaire pour l'inference de rafale m_axi.
        }
        if (lastE) break;
    }
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

static void thr_count_sink(
    hls::stream<ap_uint<64>>& total_i,
    ap_uint<64>*              nMinizrs)
{
#pragma HLS INLINE off
    *nMinizrs = total_i.read();
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
void minimizer(
    const ap_uint<MEM_WIDTH>*    packed_sequence,
    ap_uint<GROUP_BITS>*         tab_hash,
    ap_uint<64>*                 nMinizrs,
    ap_uint<64>                  n_bases
){

#pragma HLS INTERFACE m_axi port=packed_sequence offset=slave bundle=gmem0 num_read_outstanding=32  depth=COSIM_MAX_BEATS max_read_burst_length=64
#pragma HLS INTERFACE m_axi port=tab_hash        offset=slave bundle=gmem1 num_write_outstanding=32 depth=COSIM_MAX_HASHGROUP max_write_burst_length=128
#pragma HLS INTERFACE s_axilite port=nMinizrs
#pragma HLS INTERFACE s_axilite port=n_bases
#pragma HLS INTERFACE s_axilite port=return

#pragma HLS DATAFLOW
    constexpr int W8 = 8 * PAR_FACTOR;
    constexpr int W2 = 2 * PAR_FACTOR;
    constexpr int WS = 2 * SMER_SIZE * PAR_FACTOR;
    hls::stream<ap_uint<MEM_WIDTH>> s_mem;
    hls::stream<ap_uint<W8>> s_base_raw;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_raw_v;
    hls::stream<ap_uint<W2>> s_base_2b;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_2b_v;
    hls::stream<ap_uint<W2>> s_base_al;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_al_v;
    hls::stream<ap_uint<WS>> s_smer;
    hls::stream<ap_uint<PAR_FACTOR>> s_smer_v;
    hls::stream<ap_uint<WS>> s_hash;
    hls::stream<ap_uint<PAR_FACTOR>> s_hash_v;
    hls::stream<ap_uint<WS>> s_win;
    hls::stream<ap_uint<PAR_FACTOR>> s_win_v;
    hls::stream<ap_uint<WS>> s_min;
    hls::stream<ap_uint<PAR_FACTOR>> s_min_v;
    hls::stream<ap_uint<WS>> s_ded;
    hls::stream<ap_uint<PAR_FACTOR>> s_ded_v;
    hls::stream<ap_uint<64*PAR_FACTOR>> s_elem;
    hls::stream<ap_uint<8>>             s_elem_count;
    hls::stream<ap_uint<GROUP_BITS>>    s_flat_v;
    hls::stream<bool>                  s_flat_e;
    hls::stream<ap_uint<64>>            s_total;

#pragma HLS STREAM variable=s_mem        depth=2
#pragma HLS STREAM variable=s_base_raw   depth=2
#pragma HLS STREAM variable=s_base_raw_v depth=2
#pragma HLS STREAM variable=s_base_2b    depth=2
#pragma HLS STREAM variable=s_base_2b_v  depth=2
#pragma HLS STREAM variable=s_base_al    depth=2
#pragma HLS STREAM variable=s_base_al_v  depth=2
#pragma HLS STREAM variable=s_smer       depth=2
#pragma HLS STREAM variable=s_smer_v     depth=2
#pragma HLS STREAM variable=s_hash       depth=2
#pragma HLS STREAM variable=s_hash_v     depth=2
#pragma HLS STREAM variable=s_win        depth=2
#pragma HLS STREAM variable=s_win_v      depth=2
#pragma HLS STREAM variable=s_min        depth=2
#pragma HLS STREAM variable=s_min_v      depth=2
#pragma HLS STREAM variable=s_ded        depth=2
#pragma HLS STREAM variable=s_ded_v      depth=2
#pragma HLS STREAM variable=s_elem       depth=64
#pragma HLS STREAM variable=s_elem_count depth=64
#pragma HLS STREAM variable=s_flat_v     depth=64
#pragma HLS STREAM variable=s_flat_e     depth=64
#pragma HLS STREAM variable=s_total      depth=2

    thr_reader_mem<PAR_FACTOR>(packed_sequence, n_bases, s_mem);
    thr_unpack<PAR_FACTOR>(s_mem, n_bases, s_base_raw, s_base_raw_v);
    thr_adapter_hls<PAR_FACTOR>(s_base_raw, s_base_raw_v, s_base_2b, s_base_2b_v);
    m3_thr_base_adapter<PAR_FACTOR>(s_base_2b, s_base_2b_v, s_base_al, s_base_al_v);
    thr_smer_gen<PAR_FACTOR>(s_base_al, s_base_al_v, s_smer, s_smer_v);
    thr_hash<PAR_FACTOR>(s_smer, s_smer_v, s_hash, s_hash_v);
    thr_adapter_smer<PAR_FACTOR>(s_hash, s_hash_v, s_win, s_win_v);
    thr_min_v8<PAR_FACTOR>(s_win, s_win_v, s_min, s_min_v);
    thr_dedup_v8<PAR_FACTOR>(s_min, s_min_v, s_ded, s_ded_v);
    thr_compact<PAR_FACTOR>(s_ded, s_ded_v, s_elem, s_elem_count);
    thr_flatten<PAR_FACTOR>(s_elem, s_elem_count, s_flat_v, s_flat_e, s_total);
    thr_write_burst(s_flat_v, s_flat_e, tab_hash, n_bases);
    thr_count_sink(s_total, nMinizrs);
}
