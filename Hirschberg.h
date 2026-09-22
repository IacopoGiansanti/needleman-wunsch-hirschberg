#ifndef HIRSCHBERG_H
#define HIRSCHBERG_H

#include "alignment.h"
#include "scoring.h"
#include "strand.h"

Alignment hirschberg(const Strand *v, const Strand *w, const scoringModel *model); 

#endif // !HIRSCHBERG_H
