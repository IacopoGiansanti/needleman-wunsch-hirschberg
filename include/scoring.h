#ifndef SCORING_H
#define SCORING_H

#include "base.h"

typedef struct {
    int match;
    int mismatch;
    int indel;
} scoringModel;

int scoring(Base a, Base b, const scoringModel *model);

#endif // !SCORING_H
