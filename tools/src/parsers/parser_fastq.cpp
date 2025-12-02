#include "parser_fastq.hpp"
#include <cstdio>
#include <cstdlib>

#define BUFFER_SIZE 8192

bool parser_fastq::execute(stream_reader&  i_file,  stream_writer&  o_file)
{
    char buffer[BUFFER_SIZE];
    size_t n;
    int state = 1;  // 1 = header, 2 = sequence, 3 = '+', 4 = quality
    int inside_line = 0;

    // Pour reconstituer la ligne de séquence
    char *seq = NULL;
    size_t seq_len = 0, seq_cap = 0;
    while ((n = fread(buffer, 1, BUFFER_SIZE, f)) > 0) {
        for (size_t i = 0; i < n; i++) {
            char c = buffer[i];

            if (!inside_line) {  
                inside_line = 1;
                if (state == 2) { 
                    seq_len = 0;
                }
            }

            if (c == '\n' || c == '\r') {
                inside_line = 0;
                if (state == 2) {
                    // Fin de ligne de séquence → on l'affiche
                    if (seq_len > 0) {
                        seq[seq_len] = '\0';
                        printf("%s\n", seq);
                    }
                }
                state = (state % 4) + 1;
                continue;
            }

            // Enregistrer la ligne de séquence uniquement
            if (state == 2) {
                if (seq_len + 1 >= seq_cap) {
                    seq_cap = seq_cap ? seq_cap * 2 : 256;
                    seq = realloc(seq, seq_cap);
                }
                seq[seq_len++] = c;
            }
        }
    }
    free(seq);

    return true;
}
