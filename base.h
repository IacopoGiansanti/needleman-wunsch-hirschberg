#ifndef BASE_H
#define BASE_H

#include <stdint.h>

typedef uint8_t Base;

enum {
    BASE_A = 0,
    BASE_T = 1,
    BASE_G = 2,
    BASE_C = 3,
    BASE_GAP = 4
};

Base  char_to_base(char c);
char  base_to_char(Base b);

#endif
