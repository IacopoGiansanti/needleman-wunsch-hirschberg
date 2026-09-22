#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "alignment.h"
#include "base.h"
#include "fasta.h"
#include "Hirschberg.h"
#include "direction.h"
#include "Needleman-Wunsch.h"
#include "scoring.h"
#include "strand.h"

#define LINE_WIDTH 80
#define DEFAULT_SEQ1 "fasta/HIV-1-95TNIH047.fasta"
#define DEFAULT_SEQ2 "fasta/HIV-1-95TNIH022.fasta"
#define DEFAULT_REPORT "results/comparison.txt"

typedef struct {
    size_t matches;
    size_t mismatches;
    size_t gaps;
} AlignmentStats;

typedef struct {
    size_t working_bytes;   /* memoria temporanea dell'algoritmo */
    size_t output_bytes;    /* memoria necessaria per l'Alignment restituito */
    size_t peak_bytes;      /* working + output */
} MemoryStats;

static double bytes_to_mib(size_t bytes) {
    return (double)bytes / (1024.0 * 1024.0);
}

/*
 * Stima esatta del payload heap allocato dalla nostra implementazione,
 * esclusi:
 *   - input v e w, comuni ai due algoritmi;
 *   - overhead interno di malloc/libc;
 *   - stack.
 *
 * Needleman-Wunsch mantiene contemporaneamente:
 *   score[(n+1) x (m+1)] di int
 *   pred [(n+1) x (m+1)] di Direction
 * più le due sequenze dell'allineamento risultante.
 */
static MemoryStats nw_memory_stats(size_t n, size_t m) {
    size_t cells = (n + 1) * (m + 1);

    size_t working =
        cells * sizeof(int) +
        cells * sizeof(Direction);

    size_t output =
        2 * sizeof(Strand) +
        2 * (n + m) * sizeof(Base);

    MemoryStats stats = {
        .working_bytes = working,
        .output_bytes = output,
        .peak_bytes = working + output
    };

    return stats;
}

/*
 * Nella nostra Hirschberg, al picco della chiamata di livello superiore,
 * mentre viene calcolato suffix_next, sono vive contemporaneamente
 * quattro colonne di (n+1) int:
 *
 *   prefix
 *   suffix
 *   previous di nw_score_reverse()
 *   current  di nw_score_reverse()
 *
 * Gli array temporanei vengono liberati prima delle chiamate ricorsive,
 * quindi non si sommano lungo lo stack ricorsivo.
 */
static MemoryStats hb_memory_stats(size_t n, size_t m) {
    (void)m;

    size_t working =
        4 * (n + 1) * sizeof(int);

    size_t output =
        2 * sizeof(Strand) +
        2 * (n + m) * sizeof(Base);

    MemoryStats stats = {
        .working_bytes = working,
        .output_bytes = output,
        .peak_bytes = working + output
    };

    return stats;
}

static double elapsed_seconds(clock_t start, clock_t end) {
    return (double)(end - start) / (double)CLOCKS_PER_SEC;
}

static AlignmentStats alignment_stats(const Alignment *alignment) {
    AlignmentStats stats = {0, 0, 0};
    size_t length = alignment->sequence1->length;

    for(size_t i = 0; i < length; i++) {
        Base a = alignment->sequence1->bases[i];
        Base b = alignment->sequence2->bases[i];

        if(a == BASE_GAP || b == BASE_GAP)
            stats.gaps++;
        else if(a == b)
            stats.matches++;
        else
            stats.mismatches++;
    }

    return stats;
}

static bool alignment_is_valid(const Alignment *alignment,
                               const Strand *v,
                               const Strand *w,
                               const scoringModel *model) {
    if(alignment == NULL || alignment->sequence1 == NULL || alignment->sequence2 == NULL)
        return false;

    if(alignment->sequence1->length != alignment->sequence2->length)
        return false;

    size_t i = 0;
    size_t j = 0;
    int score = 0;

    for(size_t k = 0; k < alignment->sequence1->length; k++) {
        Base a = alignment->sequence1->bases[k];
        Base b = alignment->sequence2->bases[k];

        if(a == BASE_GAP && b == BASE_GAP)
            return false;

        if(a != BASE_GAP) {
            if(i >= v->length || a != v->bases[i])
                return false;
            i++;
        }

        if(b != BASE_GAP) {
            if(j >= w->length || b != w->bases[j])
                return false;
            j++;
        }

        score += scoring(a, b, model);
    }

    return i == v->length &&
           j == w->length &&
           score == alignment->score;
}

static bool same_alignment(const Alignment *a, const Alignment *b) {
    if(a->sequence1->length != b->sequence1->length ||
       a->sequence2->length != b->sequence2->length)
        return false;

    size_t length = a->sequence1->length;

    return memcmp(a->sequence1->bases, b->sequence1->bases,
                  length * sizeof(Base)) == 0 &&
           memcmp(a->sequence2->bases, b->sequence2->bases,
                  length * sizeof(Base)) == 0;
}

static void write_alignment_block(FILE *file,
                                  const char *title,
                                  const Alignment *alignment) {
    fprintf(file, "\n============================================================\n");
    fprintf(file, "%s\n", title);
    fprintf(file, "============================================================\n\n");

    size_t length = alignment->sequence1->length;

    for(size_t start = 0; start < length; start += LINE_WIDTH) {
        size_t end = start + LINE_WIDTH;
        if(end > length)
            end = length;

        fprintf(file, "[%zu-%zu]\n", start + 1, end);

        fprintf(file, "seq1: ");
        for(size_t i = start; i < end; i++)
            fputc(base_to_char(alignment->sequence1->bases[i]), file);
        fputc('\n', file);

        fprintf(file, "      ");
        for(size_t i = start; i < end; i++) {
            Base a = alignment->sequence1->bases[i];
            Base b = alignment->sequence2->bases[i];

            if(a == BASE_GAP || b == BASE_GAP)
                fputc(' ', file);
            else if(a == b)
                fputc('|', file);
            else
                fputc('.', file);
        }
        fputc('\n', file);

        fprintf(file, "seq2: ");
        for(size_t i = start; i < end; i++)
            fputc(base_to_char(alignment->sequence2->bases[i]), file);

        fprintf(file, "\n\n");
    }
}

static int write_report(const char *filename,
                        const char *seq1_path,
                        const char *seq2_path,
                        const Strand *v,
                        const Strand *w,
                        const scoringModel *model,
                        const Alignment *nw,
                        double nw_time,
                        MemoryStats nw_memory,
                        const Alignment *hb,
                        double hb_time,
                        MemoryStats hb_memory,
                        bool nw_valid,
                        bool hb_valid) {
    FILE *file = fopen(filename, "w");
    if(file == NULL) {
        perror(filename);
        return -1;
    }

    AlignmentStats nw_stats = alignment_stats(nw);
    AlignmentStats hb_stats = alignment_stats(hb);

    bool same_score = nw->score == hb->score;
    bool same_path = same_alignment(nw, hb);

    fprintf(file, "GLOBAL ALIGNMENT COMPARISON\n");
    fprintf(file, "============================================================\n\n");

    fprintf(file, "INPUT\n");
    fprintf(file, "  sequence 1 : %s\n", seq1_path);
    fprintf(file, "  length     : %zu nt\n", v->length);
    fprintf(file, "  sequence 2 : %s\n", seq2_path);
    fprintf(file, "  length     : %zu nt\n\n", w->length);

    fprintf(file, "SCORING MODEL\n");
    fprintf(file, "  match      : %+d\n", model->match);
    fprintf(file, "  mismatch   : %+d\n", model->mismatch);
    fprintf(file, "  gap        : %+d\n\n", model->indel);

    fprintf(file, "SUMMARY\n");
    fprintf(file, "%-20s %10s %10s %10s %10s %10s %12s\n",
            "algorithm", "score", "length", "matches", "mismatch", "gaps", "CPU time(s)");
    fprintf(file, "%-20s %10d %10zu %10zu %10zu %10zu %12.6f\n",
            "Needleman-Wunsch", nw->score, nw->sequence1->length,
            nw_stats.matches, nw_stats.mismatches, nw_stats.gaps, nw_time);
    fprintf(file, "%-20s %10d %10zu %10zu %10zu %10zu %12.6f\n",
            "Hirschberg", hb->score, hb->sequence1->length,
            hb_stats.matches, hb_stats.mismatches, hb_stats.gaps, hb_time);

    fprintf(file, "\nMEMORY (peak heap payload, algorithm-specific)\n");
    fprintf(file, "%-20s %14s %14s %14s\n",
            "algorithm", "working MiB", "output MiB", "peak MiB");
    fprintf(file, "%-20s %14.3f %14.3f %14.3f\n",
            "Needleman-Wunsch",
            bytes_to_mib(nw_memory.working_bytes),
            bytes_to_mib(nw_memory.output_bytes),
            bytes_to_mib(nw_memory.peak_bytes));
    fprintf(file, "%-20s %14.3f %14.3f %14.3f\n",
            "Hirschberg",
            bytes_to_mib(hb_memory.working_bytes),
            bytes_to_mib(hb_memory.output_bytes),
            bytes_to_mib(hb_memory.peak_bytes));

    fprintf(file, "\n  NW/Hirschberg peak ratio : %.2fx\n",
            hb_memory.peak_bytes == 0
                ? 0.0
                : (double)nw_memory.peak_bytes / (double)hb_memory.peak_bytes);

    fprintf(file,
            "  NOTE: values above count heap payload allocated explicitly by the\n"
            "        current implementations. Input strands, malloc metadata and\n"
            "        stack usage are excluded.\n");

    fprintf(file, "\nCHECKS\n");
    fprintf(file, "  Needleman-Wunsch valid : %s\n", nw_valid ? "YES" : "NO");
    fprintf(file, "  Hirschberg valid       : %s\n", hb_valid ? "YES" : "NO");
    fprintf(file, "  same optimal score     : %s\n", same_score ? "YES" : "NO");
    fprintf(file, "  same exact alignment   : %s\n", same_path ? "YES" : "NO");

    if(same_score && !same_path) {
        fprintf(file,
                "\nNOTE: the two algorithms have the same optimal score but produced\n"
                "different optimal alignments. This is possible when the dynamic\n"
                "programming graph contains ties.\n");
    }

    write_alignment_block(file, "NEEDLEMAN-WUNSCH ALIGNMENT", nw);
    write_alignment_block(file, "HIRSCHBERG ALIGNMENT", hb);

    if(fclose(file) != 0) {
        perror("fclose");
        return -1;
    }

    return 0;
}

static void print_usage(const char *program) {
    fprintf(stderr,
            "Usage:\n"
            "  %s\n"
            "  %s <sequence1.fasta> <sequence2.fasta> [report.txt]\n",
            program, program);
}

int main(int argc, char **argv) {
    const char *seq1_path = DEFAULT_SEQ1;
    const char *seq2_path = DEFAULT_SEQ2;
    const char *report_path = DEFAULT_REPORT;

    if(argc == 3 || argc == 4) {
        seq1_path = argv[1];
        seq2_path = argv[2];
        if(argc == 4)
            report_path = argv[3];
    } else if(argc != 1) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    scoringModel model = {
        .match = 1,
        .mismatch = -1,
        .indel = -1
    };

    Strand *v = fasta_read_strand(seq1_path);
    Strand *w = fasta_read_strand(seq2_path);

    clock_t start = clock();
    Alignment nw = needleman_wunsch(v, w, &model);
    clock_t end = clock();
    double nw_time = elapsed_seconds(start, end);

    start = clock();
    Alignment hb = hirschberg(v, w, &model);
    end = clock();
    double hb_time = elapsed_seconds(start, end);

    MemoryStats nw_memory = nw_memory_stats(v->length, w->length);
    MemoryStats hb_memory = hb_memory_stats(v->length, w->length);

    bool nw_valid = alignment_is_valid(&nw, v, w, &model);
    bool hb_valid = alignment_is_valid(&hb, v, w, &model);
    bool same_score = nw.score == hb.score;
    bool same_path = same_alignment(&nw, &hb);

    printf("\nGLOBAL ALIGNMENT COMPARISON\n");
    printf("============================================================\n");
    printf("sequence lengths : %zu x %zu\n", v->length, w->length);
    printf("scoring          : match=%+d mismatch=%+d gap=%+d\n\n",
           model.match, model.mismatch, model.indel);

    printf("%-20s %10s %10s %12s\n",
           "algorithm", "score", "length", "CPU time(s)");
    printf("%-20s %10d %10zu %12.6f\n",
           "Needleman-Wunsch", nw.score, nw.sequence1->length, nw_time);
    printf("%-20s %10d %10zu %12.6f\n",
           "Hirschberg", hb.score, hb.sequence1->length, hb_time);

    printf("\nmemory (peak heap payload, algorithm-specific):\n");
    printf("%-20s %14s %14s %14s\n",
           "algorithm", "working MiB", "output MiB", "peak MiB");
    printf("%-20s %14.3f %14.3f %14.3f\n",
           "Needleman-Wunsch",
           bytes_to_mib(nw_memory.working_bytes),
           bytes_to_mib(nw_memory.output_bytes),
           bytes_to_mib(nw_memory.peak_bytes));
    printf("%-20s %14.3f %14.3f %14.3f\n",
           "Hirschberg",
           bytes_to_mib(hb_memory.working_bytes),
           bytes_to_mib(hb_memory.output_bytes),
           bytes_to_mib(hb_memory.peak_bytes));
    printf("NW/Hirschberg peak ratio: %.2fx\n",
           hb_memory.peak_bytes == 0
               ? 0.0
               : (double)nw_memory.peak_bytes / (double)hb_memory.peak_bytes);

    printf("\nchecks:\n");
    printf("  NW valid             : %s\n", nw_valid ? "YES" : "NO");
    printf("  Hirschberg valid     : %s\n", hb_valid ? "YES" : "NO");
    printf("  same optimal score   : %s\n", same_score ? "YES" : "NO");
    printf("  same exact alignment : %s\n", same_path ? "YES" : "NO");

    int write_status = write_report(report_path,
                                    seq1_path, seq2_path,
                                    v, w, &model,
                                    &nw, nw_time, nw_memory,
                                    &hb, hb_time, hb_memory,
                                    nw_valid, hb_valid);

    if(write_status == 0)
        printf("\nfull report written to: %s\n", report_path);

    alignment_free(&nw);
    alignment_free(&hb);
    strand_free(v);
    strand_free(w);

    if(write_status != 0 || !nw_valid || !hb_valid || !same_score)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

