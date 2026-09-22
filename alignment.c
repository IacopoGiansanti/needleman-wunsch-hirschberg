#include <stdio.h>

#include "alignment.h"
#include "base.h"

void printAlignment(Alignment a) {
    printf("final score: %d\n", a.score);

    size_t l = a.sequence1->length;

    for(size_t i=0; i<l; i++) {
        printf("%c -- %c\n", base_to_char(a.sequence1->bases[i]),base_to_char(a.sequence2->bases[i]));
    }
}

void alignment_free(Alignment *a) {
    if(a == NULL) return;

    strand_free(a->sequence1);
    strand_free(a->sequence2);

    a->sequence1 = NULL;
    a->sequence2 = NULL;
}
