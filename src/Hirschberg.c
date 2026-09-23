#include <stdlib.h>
#include <stdio.h>

#include "alignment.h"
#include "scoring.h"
#include "middle_edge.h"
#include "Hirschberg.h"

static int max(int a, int b, int c) {
    if(a >= b) 
        return a * (a >= c) + c * (a < c);
    else 
        return b * (b >= c) + c * (b < c);
}

static int *nw_score(const Strand *v, size_t v_source, size_t v_sink, 
                     const Strand *w, size_t w_source, size_t w_sink, 
                     const scoringModel *model) {

    size_t n = v_sink - v_source; 
    size_t m = w_sink - w_source; 

    int *previous = malloc((n+1) * sizeof(int));
    int *current  = malloc((n+1) * sizeof(int));
    if(!previous || !current) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    previous[0] = 0;
    for(size_t i=1; i<=n; i++) {
        previous[i] = previous[i-1] + scoring(v->bases[v_source + i - 1],BASE_GAP,model);
    }
    current[0] = 0;

    for(size_t j=1; j<=m; j++) {
        current[0] = previous[0] + scoring(BASE_GAP,w->bases[w_source + j - 1],model);
        for(size_t i=1; i<=n; i++) {
            int up       = current[i-1] + scoring(v->bases[v_source + i - 1],BASE_GAP,model);
            int left     = previous[i] + scoring(BASE_GAP,w->bases[w_source + j - 1],model);
            int diagonal = previous[i-1] + scoring(v->bases[v_source + i - 1],w->bases[w_source + j - 1],model);
            current[i] = max(up,left,diagonal);
        }
        int *temp = previous;
        previous  = current;
        current   = temp;
    }

    free(current);
    return previous;
}

static int *nw_score_reverse(const Strand *v, size_t v_source, size_t v_sink, 
                             const Strand *w, size_t w_source, size_t w_sink, 
                             const scoringModel *model) {

    size_t n = v_sink - v_source;
    size_t m = w_sink - w_source;

    int *previous = malloc((n+1) * sizeof(int));
    int *current  = malloc((n+1) * sizeof(int));
    if(!previous || !current) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    previous[n] = 0;
    for(size_t i = n; i --> 0; ) {
        previous[i] = previous[i+1] + scoring(v->bases[v_source + i], BASE_GAP, model);
    }

    for(size_t j = m; j --> 0; ) {
        current[n] = previous[n] + scoring(BASE_GAP,w->bases[w_source + j],model);
        for(size_t i = n; i --> 0; ) {
            int down = current[i+1] + scoring(v->bases[v_source + i],BASE_GAP,model);
            int right = previous[i] + scoring(BASE_GAP,w->bases[w_source + j],model);
            int diagonal = previous[i+1] + scoring(v->bases[v_source + i],w->bases[w_source + j],model);
            current[i] = max(down,right,diagonal);
        }
        int *temp = previous;
        previous = current;
        current = temp;
    }

    free(current);
    return previous;

}

static void fill(Alignment *alignment, size_t *k, Base a, Base b, const scoringModel *model) {
    alignment->sequence1->bases[*k] = a;
    alignment->sequence2->bases[*k] = b;
    alignment->score += scoring(a,b,model);

    (*k)++;
} 

static void emit(Alignment *alignment, size_t *k, const MiddleEdge *mid_edge, const Strand *v, const Strand *w, const scoringModel *model) {

    switch(mid_edge->edge) {

        case EDGE_DIAGONAL:
            fill(alignment,k,v->bases[mid_edge->from_i],w->bases[mid_edge->from_j],model);
            break;

        case EDGE_VERTICAL:
            fill(alignment,k,v->bases[mid_edge->from_i],BASE_GAP,model);
            break;

        case EDGE_HORIZONTAL:
            fill(alignment,k,BASE_GAP,w->bases[mid_edge->from_j],model);
            break;
    }
}

static void hirschberg_rec(const Strand *v, size_t v_source, size_t v_sink, const Strand *w, size_t w_source, size_t w_sink, const scoringModel *model, Alignment *alignment, size_t *length) {

    if(w_source == w_sink) {
        for(size_t i=v_source; i<v_sink; i++) {
            fill(alignment,length,v->bases[i],BASE_GAP,model);
        }
        return;
    } else if(v_source == v_sink) {
        for(size_t j=w_source; j<w_sink; j++) {
            fill(alignment,length,BASE_GAP,w->bases[j],model);
        }
        return;
    }

    size_t n = v_sink - v_source;
    size_t m = w_sink - w_source;

    size_t mid_j = w_source + m/2;

    int *prefix = nw_score(v,v_source,v_sink,w,w_source,mid_j,model);
    int *suffix = nw_score_reverse(v,v_source,v_sink,w,mid_j, w_sink,model);

    size_t mid_offset = 0;
    int max = prefix[0] + suffix[0];
    for(size_t i=1; i<=n; i++) {
        int score = prefix[i] + suffix[i];
        //length[i] = prefix[i] + suffix[i]
        if(max < score) {
            max = score;
            mid_offset = i;
        }
    }
    size_t mid_i = v_source + mid_offset;
    int *suffix_next = nw_score_reverse(v,v_source,v_sink,w,mid_j + 1,w_sink,model);

    MiddleEdge edge = mid_edge(v_source,v_sink,mid_i,mid_j,suffix,suffix_next,v,w,model);  

    free(prefix);
    free(suffix);
    free(suffix_next);

    hirschberg_rec(v, v_source, mid_i, w, w_source, mid_j, model, alignment, length); // sinistra
    emit(alignment, length, &edge, v, w, model);
    hirschberg_rec(v, edge.to_i, v_sink, w, edge.to_j, w_sink, model, alignment, length);     // destra
    
}

Alignment hirschberg(const Strand *v, const Strand *w, const scoringModel *model) {
    
    size_t max_length = v->length + w->length;

    Alignment alignment =  {
        .sequence1 = strand_create(max_length),
        .sequence2 = strand_create(max_length),
        .score     = 0
    };

    size_t length = 0;

    const Strand *rows = v;
    const Strand *cols = w;
    int flag = 0;

    if(v->length > w->length) {
        rows = w;
        cols = v;
        flag = 1;
    }

    hirschberg_rec(rows, 0, rows->length, cols, 0, cols->length, model, &alignment, &length);

    alignment.sequence1->length = length;
    alignment.sequence2->length = length;

    if(flag) {
        Strand *temp = alignment.sequence1;
        alignment.sequence1 = alignment.sequence2;
        alignment.sequence2 = temp;
    }

    return alignment;

}
