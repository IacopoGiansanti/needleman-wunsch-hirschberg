#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alignment.h"
#include "base.h"
#include "fasta.h"
#include "Hirschberg.h"
#include "Needleman-Wunsch.h"
#include "scoring.h"
#include "strand.h"

static void usage(const char *prog) {
    fprintf(stderr,
            "Usage: %s <nw|hirschberg> <seq1.fasta> <seq2.fasta> [match mismatch gap]\n",
            prog);
}

static void print_sequence(const Strand *s) {
    for(size_t i = 0; i < s->length; i++)
        putchar(base_to_char(s->bases[i]));
    putchar('\n');
}

int main(int argc, char **argv) {
    if(argc != 4 && argc != 7) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    scoringModel model = {
        .match = 1,
        .mismatch = -1,
        .indel = -1
    };

    if(argc == 7) {
        model.match = atoi(argv[4]);
        model.mismatch = atoi(argv[5]);
        model.indel = atoi(argv[6]);
    }

    Strand *v = fasta_read_strand(argv[2]);
    Strand *w = fasta_read_strand(argv[3]);

    Alignment a;
    if(strcmp(argv[1], "nw") == 0) {
        a = needleman_wunsch(v, w, &model);
    } else if(strcmp(argv[1], "hirschberg") == 0) {
        a = hirschberg(v, w, &model);
    } else {
        usage(argv[0]);
        strand_free(v);
        strand_free(w);
        return EXIT_FAILURE;
    }

    printf("SCORE\t%d\n", a.score);
    printf("SEQ1\t");
    print_sequence(a.sequence1);
    printf("SEQ2\t");
    print_sequence(a.sequence2);

    alignment_free(&a);
    strand_free(v);
    strand_free(w);
    return EXIT_SUCCESS;
}
