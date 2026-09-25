#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "strand.h"
#include "scoring.h"
#include "alignment.h"

Alignment needleman_wunsch(const Strand *v, const Strand *c, const scoringModel *model);

#endif // !NEEDLEMAN_WUNSCH_H
