#ifndef DIRECTION_H
#define DIRECTION_H

#include <stdint.h>

typedef uint8_t Direction;

enum {
    NONE = 0,
    UP = 1 << 0,
    LEFT = 1 << 1,
    DIAGONAL = 1 << 2
};

#endif // !DIRECTION_H
