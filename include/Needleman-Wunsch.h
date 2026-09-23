#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "base.h"
#include "strand.h"
#include "scoring.h"
#include "alignment.h"
#include "direction.h"

Alignment needleman_wunsch(const Strand *v, const Strand *c, const scoringModel *model);

#endif // !NEEDLEMAN_WUNSCH_H
