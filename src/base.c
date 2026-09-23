#include <stdio.h>
#include <stdlib.h>

#include "base.h"

Base char_to_base(char c) {
    switch (c) {
        case 'A': case 'a': return BASE_A;
        case 'T': case 't': return BASE_T;
        case 'G': case 'g': return BASE_G;
        case 'C': case 'c': return BASE_C;
        case '-': return BASE_GAP;
        default:
            fprintf(stderr, "Unknown base: '%c'\n", c);
            exit(EXIT_FAILURE);
    }
}

char base_to_char(Base b) {
    switch (b) {
        case BASE_A: return 'A';
        case BASE_T: return 'T';
        case BASE_G: return 'G';
        case BASE_C: return 'C';
        case BASE_GAP: return '-';
        default:     return '?';
    }
}
