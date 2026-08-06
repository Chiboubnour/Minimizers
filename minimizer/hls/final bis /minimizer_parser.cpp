#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include "minimizer.hpp"


template <int resu_size>
ap_uint<resu_size> hash_u64(ap_uint<64> key)
{
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
    key = (~key + (key << 21));
    key ^= key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key ^= key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key ^= key >> 28;
    key = (key + (key << 31));
    return key;
}

template<int PAR_FACTOR>
inline ap_uint<PAR_FACTOR> valid_bits(const ap_uint<64> reste)
{
    if (reste >= PAR_FACTOR) return ~ap_uint<PAR_FACTOR>(0);
    return (ap_uint<PAR_FACTOR>(1) << reste) - 1;
}

template<int PAR_FACTOR>
void fasta_parser_pwide(
    const ap_uint<8 * PAR_FACTOR> *raw_ptr_i,
    ap_uint<64>                    n_bytes_i,
    hls::stream<ap_uint<2 * PAR_FACTOR>>& base_stream_o,
    hls::stream<ap_uint<     PAR_FACTOR>>& base_valid_o,
    hls::stream<ap_uint<     PAR_FACTOR>>& base_newseq_o
){
#pragma HLS INLINE off

    ap_uint<8> in_header         = 0;
    bool       seq_start_pending = true;

    const ap_uint<64> n_words = (n_bytes_i + PAR_FACTOR - 1) / PAR_FACTOR;
    ap_uint<64> byte_cnt = n_bytes_i;

LOOP_PARSE:
    for (ap_uint<64> w = 0; w < n_words; w++) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=1000 max=1000

        const ap_uint<8 * PAR_FACTOR> in_word   = raw_ptr_i[w];
        const ap_uint<     PAR_FACTOR> in_bounds = valid_bits<PAR_FACTOR>(byte_cnt);

        ap_uint<2 * PAR_FACTOR> out_word   = 0;
        ap_uint<     PAR_FACTOR> out_valid  = 0;
        ap_uint<     PAR_FACTOR> out_newseq = 0;

        for (int i = 0; i < PAR_FACTOR; i++) {
            const ap_uint<8> c = in_word.range(8*i + 7, 8*i);

            if (c == (ap_uint<8>)'>')       { in_header = 1; seq_start_pending = true; }
            else if (c == (ap_uint<8>)'\n') in_header = 0;

            const ap_uint<8> upper = (c >= 'a' && c <= 'z') ? (ap_uint<8>)(c - 32) : c;
            const bool is_base = (upper == 'A') || (upper == 'C') || (upper == 'G') || (upper == 'T');
            const bool valid   = (!in_header) && is_base && (bool)in_bounds[i];

            out_newseq[i] = valid && seq_start_pending;
            if (valid && seq_start_pending) seq_start_pending = false;

            out_word.range(2*i+1, 2*i) = (upper >> 1) & 3;
            out_valid[i] = valid;
        }

        if (out_valid != 0) {
            base_stream_o.write(out_word);
            base_valid_o.write(out_valid);
            base_newseq_o.write(out_newseq);
        }

        byte_cnt -= PAR_FACTOR;
    }

    base_stream_o.write(0);
    base_valid_o.write(0);
    base_newseq_o.write(0);
}


template<int PAR_FACTOR>
void thr_base_compact(
    hls::stream<ap_uint<2 * PAR_FACTOR>>& in_stream,
    hls::stream<ap_uint<     PAR_FACTOR>>& in_valid,
    hls::stream<ap_uint<     PAR_FACTOR>>& in_newseq,
    hls::stream<ap_uint<2 * PAR_FACTOR>>& out_stream,
    hls::stream<ap_uint<     PAR_FACTOR>>& out_valid,
    hls::stream<ap_uint<     PAR_FACTOR>>& out_newseq
){
#pragma HLS INLINE off

    ap_uint<2 * PAR_FACTOR> buf_d  = 0;   // buffer_cnt voies valides, voie 0 = bits[1:0]
    ap_uint<     PAR_FACTOR> buf_ns = 0;
    int buffer_cnt = 0;                   // toujours < PAR_FACTOR

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=1000 max=1000

        ap_uint<2 * PAR_FACTOR> word    = in_stream.read();
        ap_uint<     PAR_FACTOR> valid  = in_valid.read();
        ap_uint<     PAR_FACTOR> newseq = in_newseq.read();

        if (valid == 0) {
            if (buffer_cnt > 0) {
                out_stream.write(buf_d);
                out_valid.write(valid_bits<PAR_FACTOR>(buffer_cnt));
                out_newseq.write(buf_ns);
            }
            out_stream.write(0);
            out_valid.write(0);
            out_newseq.write(0);
            break;
        }


        // ap_uint<8> (pas <4>) : le compte inclusif scan[PAR_FACTOR-1] peut
        // valoir PAR_FACTOR (toutes les voies valides), qui deborde <4> des
        // que PAR_FACTOR >= 16 (constate : corruption totale de la
        // compaction a PAR_FACTOR=16, le compte silencieusement tronque a 0
        // quand tout le mot est valide). <8> couvre PAR_FACTOR jusqu'a 255,
        // coherent avec le type deja utilise en sortie de thr_compact.
        ap_uint<8> scan[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++)
            scan[i] = valid[i];
        for (int d = 1; d < PAR_FACTOR; d <<= 1) {
            ap_uint<8> nxt[PAR_FACTOR];
            for (int i = 0; i < PAR_FACTOR; i++)
                nxt[i] = (i >= d) ? (ap_uint<8>)(scan[i] + scan[i - d]) : scan[i];
            for (int i = 0; i < PAR_FACTOR; i++)
                scan[i] = nxt[i];
        }
        const ap_uint<8> k = scan[PAR_FACTOR - 1];

        ap_uint<2 * PAR_FACTOR> compact_d  = 0;
        ap_uint<     PAR_FACTOR> compact_ns = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            const ap_uint<8> pos  = scan[i] - valid[i];   // prefix exclusif
            const ap_uint<2> lane = word.range(2 * i + 1, 2 * i);
            if (valid[i]) {
                compact_d  |= (ap_uint<2 * PAR_FACTOR>(lane)     << (2 * pos));
                compact_ns |= (ap_uint<PAR_FACTOR>(newseq[i])    << pos);
            }
        }

        const int total = buffer_cnt + (int)k;

        ap_uint<4 * PAR_FACTOR> combined_d =
            ap_uint<4 * PAR_FACTOR>(buf_d) |
            (ap_uint<4 * PAR_FACTOR>(compact_d) << (2 * buffer_cnt));
        ap_uint<2 * PAR_FACTOR> combined_ns =
            buf_ns | (ap_uint<2 * PAR_FACTOR>(compact_ns) << buffer_cnt);

        if (total >= PAR_FACTOR) {
            out_stream.write((ap_uint<2 * PAR_FACTOR>)combined_d);
            out_valid.write(~ap_uint<PAR_FACTOR>(0));
            out_newseq.write((ap_uint<PAR_FACTOR>)combined_ns);

            buf_d      = (ap_uint<2 * PAR_FACTOR>)(combined_d  >> (2 * PAR_FACTOR));
            buf_ns     = (ap_uint<PAR_FACTOR>)    (combined_ns >> PAR_FACTOR);
            buffer_cnt = total - PAR_FACTOR;
        } else {
            buf_d      = (ap_uint<2 * PAR_FACTOR>)combined_d;
            buf_ns     = (ap_uint<PAR_FACTOR>)combined_ns;
            buffer_cnt = total;
        }
    }
}


template<int PAR_FACTOR>
void m3_thr_base_adapter(
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& base_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_valid,
    hls::stream<ap_uint<     PAR_FACTOR >>& base_newseq,
    hls::stream<ap_uint< 2 * PAR_FACTOR >>& out_stream,
    hls::stream<ap_uint<     PAR_FACTOR >>& out_valid,
    hls::stream<ap_uint<     PAR_FACTOR >>& out_newseq
){
#pragma HLS INLINE off
    constexpr int E            = 2;
    constexpr int SHIFT        = SMER_SIZE - 1;
    constexpr int first_rounds = SHIFT / PAR_FACTOR;
    constexpr int offset       = SHIFT % PAR_FACTOR;

    ap_uint<E * PAR_FACTOR> buffer_d  = 0;
    ap_uint<    PAR_FACTOR> buffer_v  = 0;
    ap_uint<    PAR_FACTOR> buffer_ns = 0;

    for (int r = 0; r < first_rounds; r++) {
        auto d  = base_stream.read();
        auto v  = base_valid.read();
        auto ns = base_newseq.read();
        out_stream.write(d);
        out_valid.write(v);
        out_newseq.write(ns);
    }

    {
        auto d  = base_stream.read();
        auto v  = base_valid.read();
        auto ns = base_newseq.read();

        if constexpr (offset == 0) {
            out_stream.write(d);
            out_valid.write(v);
            out_newseq.write(ns);
        } else {
            ap_uint<PAR_FACTOR> first_valid = 0;
            for (int i = 0; i < offset; i++)
                first_valid[i] = 1;

            out_stream.write(d);
            out_valid.write(first_valid);
            out_newseq.write(ns);

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                buffer_d.range(E*(i+1)-1, E*i) = d.range(E*(offset+i+1)-1, E*(offset+i));
                buffer_v[i]  = v[offset + i];
                buffer_ns[i] = ns[offset + i];
            }
        }
    }

    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        auto d  = base_stream.read();
        auto v  = base_valid.read();
        auto ns = base_newseq.read();

        if (v == 0) {
            if constexpr (offset != 0) {
                if (buffer_v != 0) {
                    out_stream.write(buffer_d);
                    out_valid.write(buffer_v);
                    out_newseq.write(buffer_ns);
                }
            }
            break;
        }

        if constexpr (offset == 0) {
            out_stream.write(d);
            out_valid.write(v);
            out_newseq.write(ns);
        } else {
            ap_uint<E * PAR_FACTOR> out_data = 0;
            ap_uint<    PAR_FACTOR> out_v    = 0;
            ap_uint<    PAR_FACTOR> out_ns   = 0;

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                out_data.range(E*(i+1)-1, E*i) = buffer_d.range(E*(i+1)-1, E*i);
                out_v[i]  = buffer_v[i];
                out_ns[i] = buffer_ns[i];
            }
            for (int i = 0; i < offset; i++) {
                out_data.range(E*(PAR_FACTOR-offset+i+1)-1, E*(PAR_FACTOR-offset+i)) = d.range(E*(i+1)-1, E*i);
                out_v[PAR_FACTOR-offset+i]  = v[i];
                out_ns[PAR_FACTOR-offset+i] = ns[i];
            }

            if (out_v != 0) {
                out_stream.write(out_data);
                out_valid.write(out_v);
                out_newseq.write(out_ns);
            }

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                buffer_d.range(E*(i+1)-1, E*i) = d.range(E*(offset+i+1)-1, E*(offset+i));
                buffer_v[i]  = v[offset + i];
                buffer_ns[i] = ns[offset + i];
            }
        }
    }

    out_stream.write(0);
    out_valid.write(0);
    out_newseq.write(0);
}

template<int PAR_FACTOR>
void thr_smer_gen(
    hls::stream<ap_uint<2 *             PAR_FACTOR>>& base_stream_i,
    hls::stream<ap_uint<                PAR_FACTOR>>& base_valid_i,
    hls::stream<ap_uint<                PAR_FACTOR>>& base_newseq_i,
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_valid_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_contam_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_firstofseq_o
){
#pragma HLS INLINE off
    constexpr int SMER_BITS    = 2 * SMER_SIZE;
    constexpr int HIST_BITS    = 2 * (SMER_SIZE - 1);
    constexpr int WIN_BITS     = HIST_BITS + 2 * PAR_FACTOR;
    constexpr int NS_HIST_BITS = SMER_SIZE - 1;
    constexpr int NS_WIN_BITS  = NS_HIST_BITS + PAR_FACTOR;
    constexpr int first_rounds = (SMER_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (SMER_SIZE - 1) % PAR_FACTOR;
    ap_uint<HIST_BITS>    memory    = 0;
    ap_uint<NS_HIST_BITS> memory_ns = 0;
    ap_uint<2 * PAR_FACTOR> in_word;
    ap_uint<    PAR_FACTOR> in_valid;
    ap_uint<    PAR_FACTOR> in_newseq;
    for (int r = 0; r < first_rounds; r++) {
        in_word   = base_stream_i.read();
        in_valid  = base_valid_i.read();
        in_newseq = base_newseq_i.read();
        for (int i = 0; i < PAR_FACTOR; i++) {
            memory = (memory << 2);
            memory(1, 0) = in_word.range(2 * i + 1, 2 * i);
            memory_ns = (memory_ns << 1);
            memory_ns[0] = in_newseq[i];
        }
    }
    if constexpr (n_last_round != 0) {
        in_word   = base_stream_i.read();
        in_valid  = base_valid_i.read();
        in_newseq = base_newseq_i.read();
        for (int i = 0; i < n_last_round; i++) {
            memory = (memory << 2);
            memory(1, 0) = in_word.range(2 * i + 1, 2 * i);
            memory_ns = (memory_ns << 1);
            memory_ns[0] = in_newseq[i];
        }
    }
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        in_word   = base_stream_i.read();
        in_valid  = base_valid_i.read();
        in_newseq = base_newseq_i.read();
        if (in_valid == 0) break;
        ap_uint<2 * PAR_FACTOR> new_rev = 0;
        for (int j = 0; j < PAR_FACTOR; j++) {
            new_rev.range(2*(PAR_FACTOR-1-j)+1, 2*(PAR_FACTOR-1-j)) =
                in_word.range(2*j+1, 2*j);
        }
        ap_uint<PAR_FACTOR> new_rev_ns = 0;
        for (int j = 0; j < PAR_FACTOR; j++)
            new_rev_ns[PAR_FACTOR-1-j] = in_newseq[j];

        const ap_uint<WIN_BITS>    win    = (memory, new_rev);
        const ap_uint<NS_WIN_BITS> win_ns = (memory_ns, new_rev_ns);
        ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_out = 0;
        ap_uint<PAR_FACTOR> contam = 0;
        ap_uint<PAR_FACTOR> firstofseq = 0;
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

            const ap_uint<SMER_SIZE> fwd_ns = win_ns.range((PAR_FACTOR-1-i)+SMER_SIZE-1, (PAR_FACTOR-1-i));
            contam[i]     = (fwd_ns.range(SMER_SIZE-2, 0) != 0);
            // bit le plus ancien (t=SMER_SIZE-1) : ce k-mer est le tout
            // premier (propre) d'une nouvelle sequence.
            firstofseq[i] = fwd_ns[SMER_SIZE-1];
        }
        smer_stream_o.write(packed_out);
        smer_valid_o.write(in_valid);
        smer_contam_o.write(contam);
        smer_firstofseq_o.write(firstofseq);
        memory    = win.range(HIST_BITS - 1, 0);
        memory_ns = win_ns.range(NS_HIST_BITS - 1, 0);
    }
    smer_stream_o.write(0);
    smer_valid_o.write(0);
    smer_contam_o.write(0);
    smer_firstofseq_o.write(0);
}

template<int PAR_FACTOR>
void thr_hash(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_contam_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_firstofseq_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  hash_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_valid_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_firstofseq_o
){
#pragma HLS INLINE off
    constexpr int S_BITS = 2 * SMER_SIZE;
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        const ap_uint<2*SMER_SIZE*PAR_FACTOR> packed_smer = smer_stream_i.read();
        const ap_uint<            PAR_FACTOR> valid       = smer_valid_i.read();
        const ap_uint<            PAR_FACTOR> contam      = smer_contam_i.read();
        const ap_uint<            PAR_FACTOR> firstofseq  = smer_firstofseq_i.read();
        if (valid == 0) {
            hash_stream_o.write(0);
            hash_valid_o.write(0);
            hash_firstofseq_o.write(0);
            break;
        }
        ap_uint<2*SMER_SIZE*PAR_FACTOR> packed_hash = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            const ap_uint<S_BITS> smer = packed_smer.range((i+1)*S_BITS-1, i*S_BITS);
            const ap_uint<S_BITS> h    = contam[i] ? ~ap_uint<S_BITS>(0)
                                                    : hash_u64<S_BITS>((ap_uint<64>)smer);
            packed_hash.range((i+1)*S_BITS-1, i*S_BITS) = h;
        }
        hash_stream_o.write(packed_hash);
        hash_valid_o.write(valid);
        hash_firstofseq_o.write(firstofseq);
    }
}

template<int PAR_FACTOR>
void thr_adapter_smer(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_firstofseq_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  smer_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_valid_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  smer_firstofseq_o
){
#pragma HLS INLINE off
    constexpr int S_BITS       = 2 * SMER_SIZE;
    constexpr int SHIFT        = WINDOW_SIZE - 1;
    constexpr int first_rounds = SHIFT / PAR_FACTOR;
    constexpr int offset       = SHIFT % PAR_FACTOR;

    ap_uint<S_BITS * PAR_FACTOR> buffer_d  = 0;
    ap_uint<PAR_FACTOR>          buffer_v  = 0;
    ap_uint<PAR_FACTOR>          buffer_fs = 0;

    for (int r = 0; r < first_rounds; r++) {
        auto d  = smer_stream_i.read();
        auto v  = smer_valid_i.read();
        auto fs = smer_firstofseq_i.read();
        smer_stream_o.write(d);
        smer_valid_o.write(v);
        smer_firstofseq_o.write(fs);
    }

    {
        auto d  = smer_stream_i.read();
        auto v  = smer_valid_i.read();
        auto fs = smer_firstofseq_i.read();

        if constexpr (offset == 0) {
            smer_stream_o.write(d);
            smer_valid_o.write(v);
            smer_firstofseq_o.write(fs);
        } else {
            ap_uint<PAR_FACTOR> first_valid = 0;
            for (int i = 0; i < offset; i++)
                first_valid[i] = 1;

            smer_stream_o.write(d);
            smer_valid_o.write(first_valid);
            smer_firstofseq_o.write(fs);

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                buffer_d.range((i+1)*S_BITS-1, i*S_BITS) =
                    d.range((offset+i+1)*S_BITS-1, (offset+i)*S_BITS);
                buffer_v[i]  = v[offset + i];
                buffer_fs[i] = fs[offset + i];
            }
        }
    }

    while (true)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        auto d  = smer_stream_i.read();
        auto v  = smer_valid_i.read();
        auto fs = smer_firstofseq_i.read();

        if (v == 0) {
            if constexpr (offset != 0) {
                if (buffer_v != 0) {
                    smer_stream_o.write(buffer_d);
                    smer_valid_o.write(buffer_v);
                    smer_firstofseq_o.write(buffer_fs);
                }
            }
            break;
        }

        if constexpr (offset == 0) {
            smer_stream_o.write(d);
            smer_valid_o.write(v);
            smer_firstofseq_o.write(fs);
        } else {
            ap_uint<S_BITS*PAR_FACTOR> out_data = 0;
            ap_uint<PAR_FACTOR>        out_v    = 0;
            ap_uint<PAR_FACTOR>        out_fs   = 0;

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                out_data.range((i+1)*S_BITS-1, i*S_BITS) = buffer_d.range((i+1)*S_BITS-1, i*S_BITS);
                out_v[i]  = buffer_v[i];
                out_fs[i] = buffer_fs[i];
            }
            for (int i = 0; i < offset; i++) {
                out_data.range((PAR_FACTOR-offset+i+1)*S_BITS-1, (PAR_FACTOR-offset+i)*S_BITS) =   d.range((i+1)*S_BITS-1, i*S_BITS);
                out_v[PAR_FACTOR-offset+i]  = v[i];
                out_fs[PAR_FACTOR-offset+i] = fs[i];
            }

            if (out_v != 0) {
                smer_stream_o.write(out_data);
                smer_valid_o.write(out_v);
                smer_firstofseq_o.write(out_fs);
            }

            for (int i = 0; i < PAR_FACTOR - offset; i++) {
                buffer_d.range((i+1)*S_BITS-1, i*S_BITS) =   d.range((offset+i+1)*S_BITS-1, (offset+i)*S_BITS);
                buffer_v[i]  = v[offset + i];
                buffer_fs[i] = fs[offset + i];
            }
        }
    }

    smer_stream_o.write(0);
    smer_valid_o.write(0);
    smer_firstofseq_o.write(0);
}


template<int PAR_FACTOR>
void thr_min_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  hash_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_valid_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_firstofseq_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_o
){
#pragma HLS INLINE off
    constexpr int S_BITS       = 2 * SMER_SIZE;
    constexpr int HIST_BITS    = (WINDOW_SIZE - 1) * S_BITS;
    constexpr int WIN_BITS     = HIST_BITS + PAR_FACTOR * S_BITS;
    constexpr int FS_HIST_BITS = WINDOW_SIZE - 1;
    constexpr int FS_WIN_BITS  = FS_HIST_BITS + PAR_FACTOR;
    constexpr int first_rounds = (WINDOW_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (WINDOW_SIZE - 1) % PAR_FACTOR;
    ap_uint<HIST_BITS>    mem    = 0;
    ap_uint<FS_HIST_BITS> mem_fs = 0;
    ap_uint<S_BITS * PAR_FACTOR> in_word;
    ap_uint<         PAR_FACTOR> in_valid;
    ap_uint<         PAR_FACTOR> in_fs;
    for (int r = 0; r < first_rounds; r++) {
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        in_fs    = hash_firstofseq_i.read();
        mem.range(S_BITS*PAR_FACTOR*(r+1)-1, S_BITS*PAR_FACTOR*r) = in_word;
        mem_fs.range(PAR_FACTOR*(r+1)-1, PAR_FACTOR*r) = in_fs;
    }
    if constexpr (n_last_round != 0) {
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        in_fs    = hash_firstofseq_i.read();
        mem.range(HIST_BITS-1, S_BITS*PAR_FACTOR*first_rounds) =
            in_word.range(n_last_round*S_BITS-1, 0);
        mem_fs.range(FS_HIST_BITS-1, PAR_FACTOR*first_rounds) =
            in_fs.range(n_last_round-1, 0);
    }
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        in_fs    = hash_firstofseq_i.read();
        if (in_valid == 0) break;
        const ap_uint<WIN_BITS>    win    = (in_word, mem);
        const ap_uint<FS_WIN_BITS> win_fs = (in_fs, mem_fs);
        ap_uint<S_BITS * PAR_FACTOR> packed_out = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            const ap_uint<S_BITS> anchor = win.range((i+1)*S_BITS-1, i*S_BITS);
            ap_uint<S_BITS> m = anchor;
            bool window_invalid = (anchor == ~ap_uint<S_BITS>(0));
            for (int w = 1; w < WINDOW_SIZE; w++) {
                const ap_uint<S_BITS> v = win.range((i+w+1)*S_BITS-1, (i+w)*S_BITS);
                if (v < m) m = v;
                if (v == ~ap_uint<S_BITS>(0)) window_invalid = true;
                if (win_fs[i+w]) window_invalid = true;
            }
            if (window_invalid) m = ~ap_uint<S_BITS>(0);
            packed_out.range((i+1)*S_BITS-1, i*S_BITS) = m;
        }
        min_stream_o.write(packed_out);
        min_valid_o.write(in_valid);
        mem    = win.range(WIN_BITS-1, PAR_FACTOR*S_BITS);
        mem_fs = win_fs.range(FS_WIN_BITS-1, PAR_FACTOR);
    }
    min_stream_o.write(0);
    min_valid_o.write(0);
}


template<int PAR_FACTOR>
void thr_dedup_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  dedup_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  dedup_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * SMER_SIZE;   // un minimizer
    const ap_uint<S_BITS> SENTINEL = ~ap_uint<S_BITS>(0);

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
            scan_has[i] = in_valid[i] && (scan_val[i] != SENTINEL);
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
            const ap_uint<S_BITS> v = in_word.range((i + 1) * S_BITS - 1, i * S_BITS);
            if (in_valid[i] && v != SENTINEL) {
                const bool            prev_has = (i == 0) ? has_last : (scan_has[i - 1] || has_last);
                const ap_uint<S_BITS> prev_val = (i == 0) ? last     : (scan_has[i - 1] ? scan_val[i - 1] : last);
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


        // ap_uint<8> (pas <4>) : voir le commentaire equivalent dans
        // thr_base_compact -- scan[PAR_FACTOR-1] deborde <4> des que
        // PAR_FACTOR >= 16 (toutes les voies valides -> compte tronque a 0).
        ap_uint<8> scan[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++) {
            scan[i] = valid[i];
        }
        for (int d = 1; d < PAR_FACTOR; d <<= 1) {
            ap_uint<8> nxt[PAR_FACTOR];
            for (int i = 0; i < PAR_FACTOR; i++) {
                nxt[i] = (i >= d) ? (ap_uint<8>)(scan[i] + scan[i - d]) : scan[i];
            }
            for (int i = 0; i < PAR_FACTOR; i++) {
                scan[i] = nxt[i];
            }
        }

        const ap_uint<8> count = scan[PAR_FACTOR - 1];


        ap_uint<64 * PAR_FACTOR> out_word = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            const ap_uint<8>  pos  = scan[i] - valid[i];   // prefix exclusif
            const ap_uint<64> lane = word.range((i + 1) * S_BITS - 1, i * S_BITS);
            if (valid[i]) {
                out_word |= (ap_uint<64 * PAR_FACTOR>(lane) << (64 * pos));
            }
        }

        out_stream.write(out_word);
        out_count.write((ap_uint<8>)count);
    }
}

template<int PAR_FACTOR>
void thr_axis_pack(
    hls::stream<ap_uint<64 * PAR_FACTOR>>& in_stream,
    hls::stream<ap_uint<8>>&               in_count,
    hls::stream<ap_uint<AXIS_BITS>>&       out_axis
){
#pragma HLS INLINE off
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<64 * PAR_FACTOR> word  = in_stream.read();
        const ap_uint<8>               count = in_count.read();

        ap_uint<AXIS_BITS> packed = 0;
        packed.range(64 * PAR_FACTOR - 1, 0)         = word;
        packed.range(AXIS_BITS - 1, 64 * PAR_FACTOR) = count;

        out_axis.write(packed);
        if (count == 0) break;
    }
}


void minimizer(
    const ap_uint<8 * PAR_FACTOR>*   raw_ptr,
    ap_uint<64>                      n_bytes,
    hls::stream<ap_uint<AXIS_BITS>>& elem_axis_o
){
#pragma HLS INTERFACE m_axi     port=raw_ptr offset=slave bundle=gmem0 num_read_outstanding=32 max_read_burst_length=64
#pragma HLS INTERFACE s_axilite port=n_bytes
#pragma HLS INTERFACE s_axilite port=return
#pragma HLS INTERFACE axis      port=elem_axis_o

#pragma HLS DATAFLOW
    constexpr int W2 = 2 * PAR_FACTOR;
    constexpr int WS = 2 * SMER_SIZE * PAR_FACTOR;

    hls::stream<ap_uint<W2>> s_base_raw;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_raw_v;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_raw_ns;
    hls::stream<ap_uint<W2>> s_base;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_v;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_ns;
    hls::stream<ap_uint<W2>> s_base_al;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_al_v;
    hls::stream<ap_uint<PAR_FACTOR>> s_base_al_ns;
    hls::stream<ap_uint<WS>> s_smer;
    hls::stream<ap_uint<PAR_FACTOR>> s_smer_v;
    hls::stream<ap_uint<PAR_FACTOR>> s_smer_contam;
    hls::stream<ap_uint<PAR_FACTOR>> s_smer_fs;
    hls::stream<ap_uint<WS>> s_hash;
    hls::stream<ap_uint<PAR_FACTOR>> s_hash_v;
    hls::stream<ap_uint<PAR_FACTOR>> s_hash_fs;
    hls::stream<ap_uint<WS>> s_win;
    hls::stream<ap_uint<PAR_FACTOR>> s_win_v;
    hls::stream<ap_uint<PAR_FACTOR>> s_win_fs;
    hls::stream<ap_uint<WS>> s_min;
    hls::stream<ap_uint<PAR_FACTOR>> s_min_v;
    hls::stream<ap_uint<WS>> s_ded;
    hls::stream<ap_uint<PAR_FACTOR>> s_ded_v;
    hls::stream<ap_uint<64*PAR_FACTOR>> s_elem;
    hls::stream<ap_uint<8>>             s_elem_count;

#pragma HLS STREAM variable=s_base_raw     depth=2
#pragma HLS STREAM variable=s_base_raw_v   depth=2
#pragma HLS STREAM variable=s_base_raw_ns  depth=2
#pragma HLS STREAM variable=s_base         depth=2
#pragma HLS STREAM variable=s_base_v       depth=2
#pragma HLS STREAM variable=s_base_ns      depth=2
#pragma HLS STREAM variable=s_base_al      depth=2
#pragma HLS STREAM variable=s_base_al_v    depth=2
#pragma HLS STREAM variable=s_base_al_ns   depth=2
#pragma HLS STREAM variable=s_smer         depth=2
#pragma HLS STREAM variable=s_smer_v       depth=2
#pragma HLS STREAM variable=s_smer_contam  depth=2
#pragma HLS STREAM variable=s_smer_fs      depth=2
#pragma HLS STREAM variable=s_hash         depth=2
#pragma HLS STREAM variable=s_hash_v       depth=2
#pragma HLS STREAM variable=s_hash_fs      depth=2
#pragma HLS STREAM variable=s_win          depth=2
#pragma HLS STREAM variable=s_win_v        depth=2
#pragma HLS STREAM variable=s_win_fs       depth=2
#pragma HLS STREAM variable=s_min          depth=2
#pragma HLS STREAM variable=s_min_v        depth=2
#pragma HLS STREAM variable=s_ded          depth=2
#pragma HLS STREAM variable=s_ded_v        depth=2
#pragma HLS STREAM variable=s_elem         depth=64
#pragma HLS STREAM variable=s_elem_count   depth=64

    fasta_parser_pwide<PAR_FACTOR>(raw_ptr, n_bytes, s_base_raw, s_base_raw_v, s_base_raw_ns);
    thr_base_compact<PAR_FACTOR>(s_base_raw, s_base_raw_v, s_base_raw_ns, s_base, s_base_v, s_base_ns);
    m3_thr_base_adapter<PAR_FACTOR>(s_base, s_base_v, s_base_ns, s_base_al, s_base_al_v, s_base_al_ns);
    thr_smer_gen<PAR_FACTOR>(s_base_al, s_base_al_v, s_base_al_ns, s_smer, s_smer_v, s_smer_contam, s_smer_fs);
    thr_hash<PAR_FACTOR>(s_smer, s_smer_v, s_smer_contam, s_smer_fs, s_hash, s_hash_v, s_hash_fs);
    thr_adapter_smer<PAR_FACTOR>(s_hash, s_hash_v, s_hash_fs, s_win, s_win_v, s_win_fs);
    thr_min_v8<PAR_FACTOR>(s_win, s_win_v, s_win_fs, s_min, s_min_v);
    thr_dedup_v8<PAR_FACTOR>(s_min, s_min_v, s_ded, s_ded_v);
    thr_compact<PAR_FACTOR>(s_ded, s_ded_v, s_elem, s_elem_count);
    thr_axis_pack<PAR_FACTOR>(s_elem, s_elem_count, elem_axis_o);
}
