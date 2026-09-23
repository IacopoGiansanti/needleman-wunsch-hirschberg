#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../alignment.h"
#include "../base.h"
#include "../scoring.h"
#include "../strand.h"
#include "bench_adapter.h"

typedef enum {
    ALG_NW,
    ALG_HIRSCHBERG
} Algorithm;

typedef enum {
    INPUT_RANDOM,
    INPUT_IDENTICAL,
    INPUT_MUTATED
} InputMode;

static double clock_seconds(clockid_t id) {
    struct timespec ts;
    if(clock_gettime(id, &ts) != 0) {
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/* PRNG deterministico: xorshift64*. Lo usiamo solo per generare input
 * riproducibili; non fa parte della regione temporizzata. */
static uint64_t rng_next(uint64_t *state) {
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * UINT64_C(2685821657736338717);
}

static Base random_base(uint64_t *state) {
    return (Base)(rng_next(state) & UINT64_C(3));
}

static void fill_random(Strand *s, uint64_t *state) {
    for(size_t i = 0; i < s->length; i++)
        s->bases[i] = random_base(state);
}

static Base mutate_base(Base original, uint64_t *state) {
    Base b;
    do {
        b = random_base(state);
    } while(b == original);
    return b;
}

static void fill_mutated(Strand *dst,
                         const Strand *src,
                         double mutation_rate,
                         uint64_t *state) {
    const uint64_t scale = UINT64_C(1000000);
    const uint64_t threshold = (uint64_t)(mutation_rate * (double)scale);

    for(size_t i = 0; i < src->length; i++) {
        uint64_t r = rng_next(state) % scale;
        dst->bases[i] = (r < threshold)
            ? mutate_base(src->bases[i], state)
            : src->bases[i];
    }
}

static size_t parse_size(const char *text, const char *name) {
    errno = 0;
    char *end = NULL;
    unsigned long long value = strtoull(text, &end, 10);

    if(errno != 0 || end == text || *end != '\0' || value == 0) {
        fprintf(stderr, "%s non valido: %s\n", name, text);
        exit(EXIT_FAILURE);
    }

    return (size_t)value;
}

static uint64_t parse_u64(const char *text, const char *name) {
    errno = 0;
    char *end = NULL;
    unsigned long long value = strtoull(text, &end, 10);

    if(errno != 0 || end == text || *end != '\0') {
        fprintf(stderr, "%s non valido: %s\n", name, text);
        exit(EXIT_FAILURE);
    }

    return (uint64_t)value;
}

static double parse_rate(const char *text) {
    errno = 0;
    char *end = NULL;
    double value = strtod(text, &end);

    if(errno != 0 || end == text || *end != '\0' || value < 0.0 || value > 1.0) {
        fprintf(stderr, "mutation_rate deve appartenere a [0, 1]\n");
        exit(EXIT_FAILURE);
    }

    return value;
}

static Algorithm parse_algorithm(const char *text) {
    if(strcmp(text, "nw") == 0)
        return ALG_NW;
    if(strcmp(text, "hirschberg") == 0 || strcmp(text, "hb") == 0)
        return ALG_HIRSCHBERG;

    fprintf(stderr, "Algoritmo non valido: %s (usa nw o hirschberg)\n", text);
    exit(EXIT_FAILURE);
}

static InputMode parse_mode(const char *text) {
    if(strcmp(text, "random") == 0)
        return INPUT_RANDOM;
    if(strcmp(text, "identical") == 0)
        return INPUT_IDENTICAL;
    if(strcmp(text, "mutated") == 0)
        return INPUT_MUTATED;

    fprintf(stderr,
            "Modalita' non valida: %s (usa random, identical o mutated)\n",
            text);
    exit(EXIT_FAILURE);
}

static const char *algorithm_name(Algorithm algorithm) {
    return algorithm == ALG_NW ? "nw" : "hirschberg";
}

static const char *mode_name(InputMode mode) {
    switch(mode) {
        case INPUT_RANDOM:    return "random";
        case INPUT_IDENTICAL: return "identical";
        case INPUT_MUTATED:   return "mutated";
    }
    return "unknown";
}

/* Controllo eseguito DOPO il timer. Verifica che l'allineamento ricostruisca
 * esattamente gli input e che lo score memorizzato sia coerente. */
static int alignment_is_valid(const Alignment *alignment,
                              const Strand *v,
                              const Strand *w,
                              const scoringModel *model) {
    if(!alignment || !alignment->sequence1 || !alignment->sequence2)
        return 0;

    if(alignment->sequence1->length != alignment->sequence2->length)
        return 0;

    size_t i = 0;
    size_t j = 0;
    int score = 0;

    for(size_t k = 0; k < alignment->sequence1->length; k++) {
        Base a = alignment->sequence1->bases[k];
        Base b = alignment->sequence2->bases[k];

        if(a == BASE_GAP && b == BASE_GAP)
            return 0;

        if(a != BASE_GAP) {
            if(i >= v->length || a != v->bases[i])
                return 0;
            i++;
        }

        if(b != BASE_GAP) {
            if(j >= w->length || b != w->bases[j])
                return 0;
            j++;
        }

        score += scoring(a, b, model);
    }

    return i == v->length &&
           j == w->length &&
           score == alignment->score;
}

static void usage(const char *program) {
    fprintf(stderr,
            "Uso:\n"
            "  %s ALGORITHM N M SEED MODE [MUTATION_RATE]\n\n"
            "ALGORITHM:\n"
            "  nw | hirschberg\n\n"
            "MODE:\n"
            "  random     due sequenze indipendenti\n"
            "  identical  sequenze identiche (richiede N=M)\n"
            "  mutated    seconda sequenza derivata dalla prima (richiede N=M)\n\n"
            "Esempi:\n"
            "  %s nw 5000 5000 42 random\n"
            "  %s hirschberg 5000 5000 42 mutated 0.05\n",
            program, program, program);
}

int main(int argc, char **argv) {
    if(argc < 6 || argc > 7) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    Algorithm algorithm = parse_algorithm(argv[1]);
    size_t n = parse_size(argv[2], "N");
    size_t m = parse_size(argv[3], "M");
    uint64_t seed = parse_u64(argv[4], "SEED");
    InputMode mode = parse_mode(argv[5]);
    double mutation_rate = argc == 7 ? parse_rate(argv[6]) : 0.05;

    if((mode == INPUT_IDENTICAL || mode == INPUT_MUTATED) && n != m) {
        fprintf(stderr, "La modalita' %s richiede N=M\n", mode_name(mode));
        return EXIT_FAILURE;
    }

    uint64_t state = seed ? seed : UINT64_C(0x9e3779b97f4a7c15);

    Strand *v = strand_create(n);
    Strand *w = strand_create(m);

    fill_random(v, &state);

    switch(mode) {
        case INPUT_RANDOM:
            fill_random(w, &state);
            break;
        case INPUT_IDENTICAL:
            memcpy(w->bases, v->bases, n * sizeof(Base));
            break;
        case INPUT_MUTATED:
            fill_mutated(w, v, mutation_rate, &state);
            break;
    }

    const scoringModel model = {
        .match = 1,
        .mismatch = -1,
        .indel = -1
    };

    /* Tutto cio' che riguarda generazione input e validazione resta fuori
     * dalla regione temporizzata. */
    double wall_start = clock_seconds(CLOCK_MONOTONIC);
    double cpu_start = clock_seconds(CLOCK_PROCESS_CPUTIME_ID);

    Alignment alignment;
    if(algorithm == ALG_NW)
        alignment = bench_run_needleman_wunsch(v, w, &model);
    else
        alignment = bench_run_hirschberg(v, w, &model);

    double cpu_end = clock_seconds(CLOCK_PROCESS_CPUTIME_ID);
    double wall_end = clock_seconds(CLOCK_MONOTONIC);

    int valid = alignment_is_valid(&alignment, v, w, &model);

    /* Nessun header: il runner lo aggiunge una sola volta al CSV finale. */
    printf("%s,%zu,%zu,%" PRIu64 ",%s,%.6f,%.9f,%.9f,%d,%d\n",
           algorithm_name(algorithm),
           n,
           m,
           seed,
           mode_name(mode),
           mutation_rate,
           wall_end - wall_start,
           cpu_end - cpu_start,
           alignment.score,
           valid);

    alignment_free(&alignment);
    strand_free(v);
    strand_free(w);

    return valid ? EXIT_SUCCESS : 2;
}

