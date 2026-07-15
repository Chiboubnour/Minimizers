#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
#include <cstdio>
#include "functions.hpp"

constexpr int MEM_WIDTH   = 64;
constexpr int PAR_FACTOR  = 8;
constexpr int SMER_SIZE   = 19;
constexpr int WINDOW_SIZE = 16;
constexpr int SMER_BITS   = 2 * SMER_SIZE;


template<int PAR_FACTOR>
inline ap_uint<PAR_FACTOR> valid_bits(const ap_uint<64> reste)
{
    if (reste >= PAR_FACTOR) return ~ap_uint<PAR_FACTOR>(0);
    return (ap_uint<PAR_FACTOR>(1) << reste) - 1;
}

template<int PAR_FACTOR>
void thr_reader(
    const ap_uint< 8 * PAR_FACTOR >*      base_ptr_i,
    const ap_uint<            64   >      n_bases_i,
    hls::stream< ap_uint< 8 * PAR_FACTOR > >& base_stream_o,
    hls::stream< ap_uint<     PAR_FACTOR > >& base_valid_o
){
#pragma HLS INLINE off

    ap_uint<64> n_words  = (n_bases_i + PAR_FACTOR - 1) / PAR_FACTOR;
    ap_uint<64> base_cnt =  n_bases_i;
    for (ap_uint<64> w = 0; w < n_words; w++)
    {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint< 8 * PAR_FACTOR > word  = base_ptr_i[w];
        const ap_uint<     PAR_FACTOR > valid = valid_bits<PAR_FACTOR>( base_cnt );

        base_stream_o.write(word);
        base_valid_o.write (valid);

        base_cnt -= PAR_FACTOR;
    }

    base_stream_o.write(0);
    base_valid_o.write (0);
}
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
        base_valid_o.write (in_valid);

        if (in_valid == 0)
        {
            break;
        }
    }
}
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

    ap_uint<E * buff_values> buffer_d = 0;
    ap_uint<    buff_values> buffer_v = 0;

    ap_uint<E * PAR_FACTOR> in_word;
    ap_uint<    PAR_FACTOR> in_valid;


    for (int r = 0; r < first_rounds; r++) {
        in_word  = base_stream.read();
        in_valid = base_valid.read();
        out_stream.write(in_word);
        out_valid.write (in_valid);
    }

    in_word  = base_stream.read();
    in_valid = base_valid.read();
    out_stream.write(in_word);
    out_valid.write (last_valid);
    buffer_d = in_word.range (E * PAR_FACTOR - 1, E * n_last_round);
    buffer_v = in_valid.range(    PAR_FACTOR - 1,     n_last_round);


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
template<int PAR_FACTOR>
void thr_smer_gen(
    hls::stream<ap_uint<2 *             PAR_FACTOR>>& base_stream_i,
    hls::stream<ap_uint<                PAR_FACTOR>>& base_valid_i,
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_o,
    hls::stream<ap_uint<                PAR_FACTOR>>& smer_valid_o
){
#pragma HLS INLINE off

    constexpr int SMER_BITS    = 2 * SMER_SIZE;                          // 2*19 : un s-mer
    constexpr int HIST_BITS    = 2 * (SMER_SIZE - 1);                    // warm-up (bits)
    constexpr int WIN_BITS     = HIST_BITS + 2 * PAR_FACTOR;             // fenetre combinee
    constexpr int first_rounds = (SMER_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (SMER_SIZE - 1) % PAR_FACTOR;

    ap_uint<HIST_BITS> memory = 0;          // SMER_SIZE-1 bases, plus recente en bas

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
    in_word  = base_stream_i.read();
    in_valid = base_valid_i.read();
    for (int i = 0; i < n_last_round; i++) {
        memory = (memory << 2);
        memory(1, 0) = in_word.range(2 * i + 1, 2 * i);
    }

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        in_word  = base_stream_i.read();
        in_valid = base_valid_i.read();

        if (in_valid == 0) {
            break;
        }

        ap_uint<2 * PAR_FACTOR> new_rev = 0;
        for (int j = 0; j < PAR_FACTOR; j++) {
#pragma HLS UNROLL
            new_rev.range(2 * (PAR_FACTOR - 1 - j) + 1, 2 * (PAR_FACTOR - 1 - j)) = in_word.range(2 * j + 1, 2 * j);
        }

        const ap_uint<WIN_BITS> win = (memory, new_rev);

        ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_out = 0;
        for (int i = 0; i < PAR_FACTOR; i++)
        {
#pragma HLS UNROLL
            const ap_uint<SMER_BITS> fwd = win.range(2 * (PAR_FACTOR - 1 - i) + 2 * SMER_SIZE - 1,
                                                     2 * (PAR_FACTOR - 1 - i));
            ap_uint<SMER_BITS> rev = 0;
            for (int t = 0; t < SMER_SIZE; t++)
            {
#pragma HLS UNROLL
                const ap_uint<2> base = fwd.range(2 * t + 1, 2 * t);
                rev.range(2 * (SMER_SIZE - 1 - t) + 1, 2 * (SMER_SIZE - 1 - t)) = base ^ ap_uint<2>(0x2);
            }
            const ap_uint<SMER_BITS> canon = (fwd < rev) ? fwd : rev;
            packed_out.range((i + 1) * SMER_BITS - 1, i * SMER_BITS) = canon;
        }

        smer_stream_o.write(packed_out);
        smer_valid_o.write (in_valid);

        memory = win.range(HIST_BITS - 1, 0);
    }

    smer_stream_o.write(0);
    smer_valid_o.write (0);
}
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

        const ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_smer = smer_stream_i.read();
        const ap_uint<                PAR_FACTOR> valid       = smer_valid_i.read();

        if (valid == 0) {
            hash_stream_o.write(0);
            hash_valid_o.write (0);
            break;
        }

        ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed_hash = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
#pragma HLS UNROLL
            const ap_uint<S_BITS> smer = packed_smer.range((i + 1) * S_BITS - 1, i * S_BITS);
            const ap_uint<S_BITS> h    = hash_u64<S_BITS>((ap_uint<64>)smer);
            packed_hash.range((i + 1) * S_BITS - 1, i * S_BITS) = h;
        }

        hash_stream_o.write(packed_hash);
        hash_valid_o.write (valid);
    }
}
template<int PAR_FACTOR>
void thr_adapter_smer(
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_i,
    hls::stream<ap_uint<PAR_FACTOR>>& smer_valid_i,
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& smer_stream_o,
    hls::stream<ap_uint<PAR_FACTOR>>& smer_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS       = 2 * SMER_SIZE;
    constexpr int SHIFT        = WINDOW_SIZE - 1;
    constexpr int first_rounds = SHIFT / PAR_FACTOR;
    constexpr int offset       = SHIFT % PAR_FACTOR;

    ap_uint<S_BITS * PAR_FACTOR> buffer_data = 0;
    ap_uint<PAR_FACTOR>          buffer_valid = 0;

    // Passage direct des blocs complets
    for(int r=0; r<first_rounds; r++){
        auto d = smer_stream_i.read();
        auto v = smer_valid_i.read();
        smer_stream_o.write(d);
        smer_valid_o.write(v);
    }

    // Initialisation du buffer
    {
        auto d = smer_stream_i.read();
        auto v = smer_valid_i.read();

        if constexpr(offset==0){
            smer_stream_o.write(d);
            smer_valid_o.write(v);
        }else{

            ap_uint<PAR_FACTOR> first_valid = 0;

            for(int i=0;i<offset;i++)
                first_valid[i]=1;

            smer_stream_o.write(d);
            smer_valid_o.write(first_valid);

            for(int i=0;i<PAR_FACTOR-offset;i++){
#pragma HLS UNROLL
                buffer_data.range((i+1)*S_BITS-1,i*S_BITS)=
                    d.range((offset+i+1)*S_BITS-1,(offset+i)*S_BITS);

                buffer_valid[i]=v[offset+i];
            }
        }
    }

    while(true){
#pragma HLS PIPELINE II=1

        auto d = smer_stream_i.read();
        auto v = smer_valid_i.read();

        if(v==0) break;

        if constexpr(offset==0){

            smer_stream_o.write(d);
            smer_valid_o.write(v);

        }else{

            ap_uint<S_BITS*PAR_FACTOR> out_data=0;
            ap_uint<PAR_FACTOR> out_valid=0;

            for(int i=0;i<PAR_FACTOR-offset;i++){
#pragma HLS UNROLL
                out_data.range((i+1)*S_BITS-1,i*S_BITS)=
                    buffer_data.range((i+1)*S_BITS-1,i*S_BITS);

                out_valid[i]=buffer_valid[i];
            }

            for(int i=0;i<offset;i++){
#pragma HLS UNROLL
                out_data.range((PAR_FACTOR-offset+i+1)*S_BITS-1,
                               (PAR_FACTOR-offset+i)*S_BITS)=
                    d.range((i+1)*S_BITS-1,i*S_BITS);

                out_valid[PAR_FACTOR-offset+i]=v[i];
            }

            if(out_valid!=0){
                smer_stream_o.write(out_data);
                smer_valid_o.write(out_valid);
            }

            for(int i=0;i<PAR_FACTOR-offset;i++){
#pragma HLS UNROLL
                buffer_data.range((i+1)*S_BITS-1,i*S_BITS)=
                    d.range((offset+i+1)*S_BITS-1,(offset+i)*S_BITS);

                buffer_valid[i]=v[offset+i];
            }
        }
    }

    smer_stream_o.write(0);
    smer_valid_o.write(0);
}
template<int PAR_FACTOR>
void thr_min_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  hash_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  hash_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_o,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_o
){
#pragma HLS INLINE off

    constexpr int S_BITS       = 2 * SMER_SIZE;                          // un hash
    constexpr int HIST_BITS    = (WINDOW_SIZE - 1) * S_BITS;             // warm-up (bits)
    constexpr int WIN_BITS     = HIST_BITS + PAR_FACTOR * S_BITS;        // fenetre combinee
    constexpr int first_rounds = (WINDOW_SIZE - 1) / PAR_FACTOR;
    constexpr int n_last_round = (WINDOW_SIZE - 1) % PAR_FACTOR;

    ap_uint<HIST_BITS> mem = 0;             // WINDOW_SIZE-1 hashes, plus ancien en bas

    ap_uint<S_BITS * PAR_FACTOR> in_word;
    ap_uint<         PAR_FACTOR> in_valid;

    //
    // Prologue : first_rounds paquets pleins + 1 paquet de n_last_round hashes.
    //
    for (int r = 0; r < first_rounds; r++) {
        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();
        mem.range(S_BITS * PAR_FACTOR * (r + 1) - 1, S_BITS * PAR_FACTOR * r) = in_word;
    }
    in_word  = hash_stream_i.read();
    in_valid = hash_valid_i.read();
    mem.range(HIST_BITS - 1, S_BITS * PAR_FACTOR * first_rounds) = in_word.range(n_last_round * S_BITS - 1, 0);

    //
    // Regime permanent.
    //
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        in_word  = hash_stream_i.read();
        in_valid = hash_valid_i.read();

        if (in_valid == 0) {
            break;
        }

        const ap_uint<WIN_BITS> win = (in_word, mem);

        ap_uint<S_BITS * PAR_FACTOR> packed_out = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
#pragma HLS UNROLL
            ap_uint<S_BITS> m = win.range((i + 1) * S_BITS - 1, i * S_BITS);
            for (int w = 1; w < WINDOW_SIZE; w++) {
                const ap_uint<S_BITS> v = win.range((i + w + 1) * S_BITS - 1, (i + w) * S_BITS);
                if (v < m) m = v;
            }
            packed_out.range((i + 1) * S_BITS - 1, i * S_BITS) = m;
        }

        min_stream_o.write(packed_out);
        min_valid_o.write (in_valid);

        mem = win.range(WIN_BITS - 1, PAR_FACTOR * S_BITS);
    }

    min_stream_o.write(0);
    min_valid_o.write (0);
}
#if 0
template<int PAR_FACTOR>
void thr_dedup_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  dedup_stream_o,
    hls::stream<ap_uint<             8             >>&   dedup_count_o
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * SMER_SIZE;

    ap_uint<S_BITS> last     = 0;
    bool            has_last = false;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<2*SMER_SIZE*PAR_FACTOR> in_word  = min_stream_i.read();
        const ap_uint<            PAR_FACTOR> in_valid = min_valid_i.read();

        if (in_valid == 0) {
            dedup_stream_o.write(0);
            dedup_count_o.write(0);
            break;
        }

        // ── Extraction — UNROLL ──────────────────────────────────────
        ap_uint<S_BITS> val[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++) {
            val[i] = in_word.range((i+1)*S_BITS - 1, i*S_BITS);
        }

        // ── Déduplication via last_pipe — UNROLL ─────────────────────
        ap_uint<S_BITS> last_pipe[PAR_FACTOR + 1];
        bool            has_pipe [PAR_FACTOR + 1];

        last_pipe[0] = last;
        has_pipe [0] = has_last;

        ap_uint<PAR_FACTOR> out_valid = 0;

        for (int i = 0; i < PAR_FACTOR; i++) {
            if (in_valid[i]) {
                bool is_new    = !has_pipe[i] || (val[i] != last_pipe[i]);
                out_valid[i]   = is_new;
                last_pipe[i+1] = is_new ? val[i] : last_pipe[i];
                has_pipe [i+1] = true;
            } else {
                last_pipe[i+1] = last_pipe[i];
                has_pipe [i+1] = has_pipe[i];
                out_valid[i]   = 0;
            }
        }

        last     = last_pipe[PAR_FACTOR];
        has_last = has_pipe [PAR_FACTOR];

        // ── Popcount — UNROLL ────────────────────────────────────────
        ap_uint<8> out_count = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (out_valid[i]) out_count++;
        }

        if (out_count == 0) continue;


        ap_uint<6> prefix[PAR_FACTOR];

        for (int i = 0; i < PAR_FACTOR; i++) {
            ap_uint<6> p = 0;
            for (int j = 0; j < i; j++) {
                if (out_valid[j]) p++;
            }
            prefix[i] = p;
        }


        ap_uint<2*SMER_SIZE*PAR_FACTOR> out_word = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (out_valid[i]) {
                for (int k = 0; k < PAR_FACTOR; k++) {
                    if (prefix[i] == k) {
                        out_word.range((k+1)*S_BITS-1, k*S_BITS) = val[i];
                    }
                }
            }
        }

        dedup_stream_o.write(out_word);
        dedup_count_o.write(out_count);
    }
}
#else
template<int PAR_FACTOR>
void thr_dedup_v8(
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  min_stream_i,
    hls::stream<ap_uint<                 PAR_FACTOR >>&  min_valid_i,
    hls::stream<ap_uint< 2 * SMER_SIZE * PAR_FACTOR >>&  dedup_stream_o,
    hls::stream<ap_uint<             8             >>&   dedup_count_o
){
#pragma HLS INLINE off

    constexpr int S_BITS = 2 * SMER_SIZE;

    ap_uint<S_BITS> last     = 0;
    bool            has_last = false;

    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<2*SMER_SIZE*PAR_FACTOR> in_word  = min_stream_i.read();
        const ap_uint<            PAR_FACTOR> in_valid = min_valid_i.read();

        if (in_valid == 0) {
            dedup_stream_o.write(0);
            dedup_count_o.write(0);
            break;
        }

        // Extraction
        ap_uint<S_BITS> val[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++)
            val[i] = in_word.range((i+1)*S_BITS-1, i*S_BITS);

        // ── Déduplication sans chaîne last_pipe ───────────────────────
        // Principe :
        //   - val[i] est doublon si égal à last (inter-paquet)
        //     OU égal à val[j] pour j < i valide (intra-paquet)
        //   - Toutes les comparaisons sont INDÉPENDANTES de last_pipe
        //   - Chemin critique = 1 comparateur 42-bit + arbre OR
        //     au lieu de 8 MUX 42-bit en chaîne

        // Étape 1 : comparaison avec last — PAR_FACTOR comparateurs parallèles
        bool eq_last[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++)
            eq_last[i] = has_last && (val[i] == last);

        // Étape 2 : comparaison intra-paquet — matrice PAR_FACTOR×PAR_FACTOR
        // Toutes indépendantes, chemin = 1 comparateur
        bool eq_prev[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++) {
            eq_prev[i] = false;
            for (int j = 0; j < i; j++) {
                if (in_valid[j] && (val[i] == val[j]))
                    eq_prev[i] = true;
            }
        }

        // Étape 3 : out_valid — arbre OR, profondeur log2
        ap_uint<PAR_FACTOR> out_valid = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (in_valid[i] && !eq_last[i] && !eq_prev[i])
                out_valid[i] = 1;
        }

        // Étape 4 : mise à jour last = dernière valeur valide du paquet
        // Réduction linéaire sur 1-bit — chemin court (pas de MUX 42-bit chaîné)
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (out_valid[i]) {
                last     = val[i];
                has_last = true;
            }
        }

        // Popcount
        ap_uint<8> out_count = 0;
        for (int i = 0; i < PAR_FACTOR; i++)
            if (out_valid[i]) out_count++;

        if (out_count == 0) continue;

        // Compactage
        ap_uint<6> prefix[PAR_FACTOR];
        for (int i = 0; i < PAR_FACTOR; i++) {
            ap_uint<6> p = 0;
            for (int j = 0; j < i; j++)
                if (out_valid[j]) p++;
            prefix[i] = p;
        }

        ap_uint<2*SMER_SIZE*PAR_FACTOR> out_word = 0;
        for (int i = 0; i < PAR_FACTOR; i++) {
            if (out_valid[i]) {
                for (int k = 0; k < PAR_FACTOR; k++) {
                    if (prefix[i] == k)
                        out_word.range((k+1)*S_BITS-1, k*S_BITS) = val[i];
                }
            }
        }

        dedup_stream_o.write(out_word);
        dedup_count_o.write(out_count);
    }
}
#endif
template<int PAR_FACTOR>
void thr_store_burst(
    hls::stream<ap_uint<2 * SMER_SIZE * PAR_FACTOR>>& minz_stream_i,
    hls::stream<ap_uint<8>>&                          minz_count_i,
    ap_uint<64>*                                     tab_hash,
    ap_uint<64>*                                     nElem)
{
#pragma HLS INLINE off

    constexpr int S_BITS     = 2 * SMER_SIZE;
    constexpr int BURST_SIZE = 64;

    ap_uint<64> cnt = 0;

    ap_uint<64> burst_buf[BURST_SIZE];
#pragma HLS ARRAY_PARTITION variable=burst_buf complete

    ap_uint<7> bcnt = 0;

STORE_LOOP:
    while (true) {
#pragma HLS PIPELINE II=1
#pragma HLS loop_tripcount min=100 max=100

        const ap_uint<2 * SMER_SIZE * PAR_FACTOR> packed = minz_stream_i.read();
        const ap_uint<8>                          count  = minz_count_i.read();

        if (count == 0)
            break;

        for (int i = 0; i < PAR_FACTOR; i++) {
            if (i < count) {
                burst_buf[bcnt] = packed.range((i + 1) * S_BITS - 1,
                                               i * S_BITS);
                bcnt++;

                if (bcnt == BURST_SIZE) {
FLUSH:
                    for (int j = 0; j < BURST_SIZE; j++) {
                        tab_hash[cnt++] = burst_buf[j];
                    }
                    bcnt = 0;
                }
            }
        }
    }

FLUSH_TAIL:
    for (ap_uint<7> j = 0; j < bcnt; j++) {
        tab_hash[cnt++] = burst_buf[j];
    }

    *nElem = cnt;
}

void minimizer(
    const ap_uint<8*PAR_FACTOR>* packed_sequence,
    ap_uint<64>*                 tab_hash,
    ap_uint<64>*                 nMinizrs,
    ap_uint<64>                  n_bases
){
#pragma HLS INTERFACE m_axi port=packed_sequence offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=tab_hash offset=slave bundle=gmem1
#pragma HLS INTERFACE s_axilite port=nMinizrs
#pragma HLS INTERFACE s_axilite port=n_bases
#pragma HLS INTERFACE s_axilite port=return
#pragma HLS DATAFLOW

    constexpr int W8 = 8 * PAR_FACTOR;
    constexpr int W2 = 2 * PAR_FACTOR;
    constexpr int WS = 2 * SMER_SIZE * PAR_FACTOR;

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
    hls::stream<ap_uint<8>>  s_ded_count;

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
#pragma HLS STREAM variable=s_ded_count  depth=2

    thr_reader<PAR_FACTOR>(packed_sequence, n_bases, s_base_raw, s_base_raw_v);
    thr_adapter_hls<PAR_FACTOR>(s_base_raw, s_base_raw_v, s_base_2b, s_base_2b_v);
    m3_thr_base_adapter<PAR_FACTOR>(s_base_2b, s_base_2b_v, s_base_al, s_base_al_v);
    thr_smer_gen<PAR_FACTOR>(s_base_al, s_base_al_v, s_smer, s_smer_v);
    thr_hash<PAR_FACTOR>(s_smer, s_smer_v, s_hash, s_hash_v);
    thr_adapter_smer<PAR_FACTOR>(s_hash, s_hash_v, s_win, s_win_v);
    thr_min_v8<PAR_FACTOR>(s_win, s_win_v, s_min, s_min_v);
    thr_dedup_v8<PAR_FACTOR>(s_min, s_min_v, s_ded, s_ded_count);
    thr_store_burst<PAR_FACTOR>(s_ded, s_ded_count, tab_hash, nMinizrs);

}


// ─────────────────────────────────────────────────────────────────────────────

extern "C"
void minimizer_v3(
    const uint64_t* seq,
    int s, int w, int n,
    uint64_t* tab_hash,
    uint64_t* nMin
){
    if (s != SMER_SIZE || w != WINDOW_SIZE){
        *nMin = 0;
        return;
    }

    // Correction : cast explicite du pointeur d'entrée
    const ap_uint<8 * PAR_FACTOR>* packed =
        reinterpret_cast<const ap_uint<8 * PAR_FACTOR>*>(seq);

    ap_uint<64>* th = reinterpret_cast<ap_uint<64>*>(tab_hash);
    ap_uint<64>  nm = 0;

    minimizer(packed, th, &nm, (ap_uint<64>)n);

    *nMin = (uint64_t)nm;
}





