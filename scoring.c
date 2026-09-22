#include "scoring.h"

int scoring(Base a, Base b, const scoringModel *model) {
   if(a == BASE_GAP || b == BASE_GAP) return model->indel;
   if(a == b) return model->match;
   return model->mismatch;
}
