#ifndef BENCH_ADAPTER_H
#define BENCH_ADAPTER_H

#include "alignment.h"
#include "scoring.h"
#include "strand.h"

/*
 * Adapter tra il benchmark e le implementazioni reali del progetto.
 *
 * Le firme rispecchiano quelle dei due algoritmi:
 * - input: due Strand in sola lettura
 * - scoring model in sola lettura
 * - output: Alignment
 */
Alignment bench_run_needleman_wunsch(
    const Strand *a,
    const Strand *b,
    const scoringModel *model
);

Alignment bench_run_hirschberg(
    const Strand *a,
    const Strand *b,
    const scoringModel *model
);

#endif
