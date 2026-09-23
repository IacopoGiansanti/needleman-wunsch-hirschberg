#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "strand.h"

Strand *strand_create(size_t length) {
    Strand *s = malloc(sizeof(Strand));
    if (!s) { perror("malloc"); exit(EXIT_FAILURE); }
    s->bases = malloc(length * sizeof(Base));
    if (!s->bases) { perror("malloc"); exit(EXIT_FAILURE); }
    s->length = length;
    return s;
}

Strand *strand_from_string(const char *seq) {
    size_t len = strlen(seq);
    Strand *s  = strand_create(len);
    for (size_t i = 0; i < len; i++)
        s->bases[i] = char_to_base(seq[i]);
    return s;
}

Strand *strand_copy(const Strand *src) {
    Strand *dst = strand_create(src->length);
    memcpy(dst->bases, src->bases, src->length * sizeof(Base));
    return dst;
}

void strand_print(const Strand *s, const char *label) {
    printf("%s (len=%zu): ", label, s->length);
    for (size_t i = 0; i < s->length; i++)
        putchar(base_to_char(s->bases[i]));
    putchar('\n');
}

void strand_free(Strand *s) {
    if (s) {
        free(s->bases);
        free(s);
    }
}
