#include <stdio.h>

#include "middle_edge.h"
#include "scoring.h"
#include "strand.h"

static int max(int a, int b, int c) {
    if(a >= b) 
        return a * (a >= c) + c * (a < c);
    else 
        return b * (b >= c) + c * (b < c);
}


MiddleEdge mid_edge(size_t top, size_t bottom, 
                    size_t mid_vertex_i, size_t mid_vertex_j,
                    const int *suffix, const int* suffix_next,
                    const Strand *v, const Strand *w, const scoringModel *model) {

    size_t i = mid_vertex_i - top;
    
    if(mid_vertex_i == bottom) {
        MiddleEdge mid_edge = {
            .from_i = mid_vertex_i,
            .from_j = mid_vertex_j,
            .to_i   = mid_vertex_i,
            .to_j   = mid_vertex_j + 1,
            .edge   = EDGE_HORIZONTAL
        };
        return mid_edge;
    }     

    int next_down     = suffix[i + 1]      + scoring(v->bases[mid_vertex_i],BASE_GAP,model);
    int next_right    = suffix_next[i]     + scoring(BASE_GAP,w->bases[mid_vertex_j],model);
    int next_diagonal = suffix_next[i + 1] + scoring(v->bases[mid_vertex_i],w->bases[mid_vertex_j],model);
    int final =  max(next_down,next_right,next_diagonal);

    MiddleEdge mid_edge = {
        .from_i = mid_vertex_i,
        .from_j = mid_vertex_j,
    };

    if(final == next_diagonal) {
        mid_edge.to_i = mid_vertex_i + 1;
        mid_edge.to_j = mid_vertex_j + 1;
        mid_edge.edge = EDGE_DIAGONAL;
    } else if(final == next_right) {
        mid_edge.to_i = mid_vertex_i;
        mid_edge.to_j = mid_vertex_j + 1;
        mid_edge.edge = EDGE_HORIZONTAL;
    } else {
        mid_edge.to_i = mid_vertex_i + 1;
        mid_edge.to_j = mid_vertex_j;
        mid_edge.edge = EDGE_VERTICAL;
    }

    return mid_edge;
}
