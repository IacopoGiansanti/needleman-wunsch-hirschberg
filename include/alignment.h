#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "strand.h"

typedef struct {
    Strand *sequence1;
    Strand *sequence2;
    int score;
} Alignment;

void printAlignment(Alignment a);
void alignment_free(Alignment *a);

#endif // !ALIGNMENT_H
