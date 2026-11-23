#include <stdlib.h>

#include "../include/direction.h"
#include "../include/strand.h"
#include "../include/base.h"

Direction graph[] = {UP, UP, DIAGONAL, LEFT, LEFT, DIAGONAL, UP, DIAGONAL, LEFT, LEFT, DIAGONAL, UP};
size_t n = 12;

typedef struct {
    Strand *v;
    Strand *w;
} Pair;

Base random_base(void) {
    return (Base)(rand() % 4);
}

Pair allineamento(void) {

    Pair result = {
        .v = strand_create(n),
        .w = strand_create(n)
    };

    size_t iv = 0;
    size_t iw = 0;
    
    for (size_t k = 0; k < n; k++) {
       switch(graph[k]) {
           case DIAGONAL: {
                              Base b = random_base();
                              result.v -> bases[iv++] = b;
                              result.w -> bases[iw++] = b;
                              break;
                          }
            case UP:
                          result.v->bases[iv++] = random_base();
                          break;
            case LEFT:
                          result.w->bases[iw++] = random_base();
                          break;
            default:
                          strand_free(result.v);
                          strand_free(result.w);
                          exit(EXIT_FAILURE);
       } 
    }

    result.v->length = iv;
    result.w->length = iw;

    strand_print(result.v, "v");
    strand_print(result.w, "w");

    return result;
}

int main() {
    Pair result = allineamento();

    strand_free(result.v);
    strand_free(result.w);

    return 0;
}
