#include "bench_adapter.h"

#include "../Needleman-Wunsch.h"
#include "../Hirschberg.h"
#include "../alignment.h"
#include "../strand.h"
#include "../scoring.h"

Alignment bench_run_needleman_wunsch(
    const Strand *a,
    const Strand *b,
    const scoringModel *model
) {
    return needleman_wunsch(a, b, model);
}

Alignment bench_run_hirschberg(
    const Strand *a,
    const Strand *b,
    const scoringModel *model
) {
    return hirschberg(a, b, model);
}

