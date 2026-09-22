#include <stdio.h>
#include <stdlib.h>

#include "base.h"
#include "strand.h"
#include "direction.h"
#include "scoring.h"
#include "alignment.h"
#include "Needleman-Wunsch.h"

static int max(int a, int b, int c) {
    if(a >= b) 
        return a * (a >= c) + c * (a < c);
    else 
        return b * (b >= c) + c * (b < c);
}

static void reverse(Base *bases, size_t length) {
    for(size_t i=0; i<length/2; i++) {
        Base temp = bases[i];
        bases[i] = bases[length - 1 - i];
        bases[length - 1 - i] = temp; 
    }
}

Alignment needleman_wunsch(const Strand *v, const Strand *w, const scoringModel *model) {
    size_t n = v->length;
    size_t m = w->length;

    Alignment alignment;

    alignment.sequence1 = strand_create(n+m);
    alignment.sequence2 = strand_create(n+m);

    Base *aligned_v = alignment.sequence1->bases;
    Base *aligned_w = alignment.sequence2->bases;

    int (*score)[m+1] = malloc((n+1) * sizeof *score);
    if(score == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    Direction (*pred)[m+1] = malloc((n+1) * sizeof *pred);
    if(pred == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    score[0][0] = 0;
    pred[0][0] = 0;

    for(size_t i=1; i<=n; i++) {
        score[i][0] = score[i-1][0] + scoring(v->bases[i-1],BASE_GAP,model);
        pred[i][0] = UP;
    }

    for(size_t j=1; j<=m; j++) {
        score[0][j] = score[0][j-1] + scoring(BASE_GAP,w->bases[j-1],model);
        pred[0][j] = LEFT;
    }

    for(size_t i=1; i<=n; i++) {
        for(size_t j=1; j<=m; j++) {
            int up = score[i-1][j] + scoring(v->bases[i-1],BASE_GAP,model);
            int left = score[i][j-1] + scoring(BASE_GAP,w->bases[j-1],model);
            int diagonal = score[i-1][j-1] + scoring(v->bases[i-1],w->bases[j-1],model);
            score[i][j] = max(up,left,diagonal);

            // salva predecessore
            pred[i][j] = NONE;
            if(score[i][j] == up) pred[i][j]|=UP;
            if(score[i][j] == left) pred[i][j]|=LEFT;
            if(score[i][j] == diagonal) pred[i][j]|=DIAGONAL;
        }
    }

    // effettua backtracking da (n,m)
    size_t length = 0;
    size_t i = n, j = m;

    while(i>0 || j>0) {
        if((pred[i][j]&DIAGONAL)==DIAGONAL) {
            aligned_v[length] = v->bases[i-1];
            aligned_w[length] = w->bases[j-1];

            length++;
            i--;
            j--;
        } else if((pred[i][j]&UP)==UP) {
            aligned_v[length] = v->bases[i-1];
            aligned_w[length] = BASE_GAP;

            length++;
            i--;
        } else {
            aligned_v[length] = BASE_GAP;
            aligned_w[length] = w->bases[j-1];

            length++;
            j--;
        }
    }

    reverse(aligned_v,length);
    reverse(aligned_w,length);

    alignment.score = score[n][m];
    alignment.sequence1->length = length;
    alignment.sequence2->length = length;

    free(score);
    free(pred);

    return alignment;
}
