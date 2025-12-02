#include "parser_fasta.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 10000
  //lit un fichier FASTA et extrait uniquement les séquences ADN, en ignorant les lignes d’en-tête (celles qui commencent par >)
 //Il reconstruit chaque séquence sur une seule ligne et l’affiche
 
bool parser_fasta::execute(stream_reader&  i_file,  stream_writer&  o_file)
{

    return true;
}



int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage : %s fichier.fasta\n", argv[0]);
        return 1;
    }

    FILE *f = fopen(argv[1], "r");
    if (!f) {
        perror("Erreur ouverture fichier");
        return 1;
    }

    char line[MAX_LINE];
    char *sequence = NULL;
    size_t seq_capacity = 0;
    size_t seq_len = 0;

    while (fgets(line, MAX_LINE, f)) {
        // En-tête FASTA
        if (line[0] == '>') {

            // Si une séquence précédente existe → on l'imprime
            if (seq_len > 0) {
                sequence[seq_len] = '\0';
                printf("%s\n", sequence);
                seq_len = 0;
            }
            continue; // On ignore la ligne d'en-tête
        }

        // Nettoyage du retour à la ligne
        line[strcspn(line, "\r\n")] = 0;

        size_t line_len = strlen(line);

        // Ajuster la taille si nécessaire
        if (seq_len + line_len + 1 > seq_capacity) {
            seq_capacity = seq_capacity * 2 + line_len + 1;
            sequence = realloc(sequence, seq_capacity);
            if (!sequence) {
                fprintf(stderr, "Erreur allocation mémoire\n");
                return 1;
            }
        }

        // Ajouter la ligne à la séquence
        memcpy(sequence + seq_len, line, line_len);
        seq_len += line_len;
    }

    // Dernière séquence du fichier
    if (seq_len > 0) {
        sequence[seq_len] = '\0';
        printf("%s\n", sequence);
    }

    free(sequence);
    fclose(f);
    return 0;
}
