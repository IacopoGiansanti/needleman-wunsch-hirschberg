#ifndef STRAND_H
#define STRAND_H

#include <stddef.h>
#include "base.h"

typedef struct {
    Base  *bases;
    size_t length;
} Strand;

Strand *strand_create(size_t length);
Strand *strand_from_string(const char *seq);
Strand *strand_copy(const Strand *src);
void    strand_print(const Strand *s, const char *label);
void    strand_free(Strand *s);

#endif
