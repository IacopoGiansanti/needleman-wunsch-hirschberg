#ifndef MIDDLE_EDGE_H
#define MIDDLE_EDGE_H

#include <stdlib.h>
#include <stdint.h>

#include "scoring.h"
#include "strand.h"

typedef uint8_t Edge; 

enum {
    EDGE_VERTICAL   = 1 << 0,
    EDGE_HORIZONTAL = 1 << 1,
    EDGE_DIAGONAL   = 1 << 2
};

typedef struct {
    size_t from_i;
    size_t from_j;
    size_t to_i;
    size_t to_j;

    Edge edge;
} MiddleEdge;

MiddleEdge mid_edge(size_t top, size_t bottom,
                    size_t mid_vertex_i, size_t mid_vertex_j,
                    const int *suffix, const int* suffix_next,
                    const Strand *v, const Strand *w, const scoringModel *model);



#endif // !MIDDLE_EDGE_H
