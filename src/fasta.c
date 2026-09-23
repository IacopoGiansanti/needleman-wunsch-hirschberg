#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "base.h"
#include "fasta.h"
#include "strand.h"

static void fasta_error(FILE *file, Base *bases, const char *filename, const char *message) {
    if (file != NULL)
        fclose(file);
    free(bases);
    fprintf(stderr, "FASTA error in '%s': %s\n", filename, message);
    exit(EXIT_FAILURE);
}

static int is_dna_base_char(int c) {
    switch (c) {
        case 'A': case 'a':
        case 'C': case 'c':
        case 'G': case 'g':
        case 'T': case 't':
            return 1;
        default:
            return 0;
    }
}

Strand *fasta_read_strand(const char *filename) {
    if (filename == NULL) {
        fprintf(stderr, "FASTA error: filename is NULL\n");
        exit(EXIT_FAILURE);
    }

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    size_t capacity = 1024;
    size_t length = 0;
    Base *bases = malloc(capacity * sizeof *bases);
    if (bases == NULL) {
        fclose(file);
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    int c;
    int seen_header = 0;
    int at_line_start = 1;

    while ((c = fgetc(file)) != EOF) {
        if (!seen_header) {
            if (at_line_start && c == '>') {
                seen_header = 1;

                /* L'intestazione non serve allo Strand: la saltiamo. */
                while ((c = fgetc(file)) != EOF && c != '\n') {
                    /* niente */
                }

                at_line_start = 1;
                continue;
            }

            if (isspace((unsigned char)c)) {
                at_line_start = (c == '\n' || c == '\r') ? 1 : at_line_start;
                continue;
            }

            fasta_error(file, bases, filename,
                        "expected a header line beginning with '>'");
        }

        if (at_line_start && c == '>') {
            fasta_error(file, bases, filename,
                        "multiple FASTA records are not supported by fasta_read_strand");
        }

        if (c == '\n' || c == '\r') {
            at_line_start = 1;
            continue;
        }

        if (isspace((unsigned char)c)) {
            /* Spazi/tabs interni alle righe vengono ignorati. */
            continue;
        }

        at_line_start = 0;

        if (!is_dna_base_char(c)) {
            char message[96];
            snprintf(message, sizeof message,
                     "unsupported sequence character '%c'", c);
            fasta_error(file, bases, filename, message);
        }

        if (length == capacity) {
            if (capacity > ((size_t)-1) / 2) {
                fasta_error(file, bases, filename, "sequence is too large");
            }

            size_t new_capacity = capacity * 2;
            Base *tmp = realloc(bases, new_capacity * sizeof *bases);
            if (tmp == NULL) {
                fclose(file);
                free(bases);
                perror("realloc");
                exit(EXIT_FAILURE);
            }

            bases = tmp;
            capacity = new_capacity;
        }

        bases[length++] = char_to_base((char)c);
    }

    if (ferror(file)) {
        fasta_error(file, bases, filename, "error while reading file");
    }

    fclose(file);

    if (!seen_header) {
        free(bases);
        fprintf(stderr, "FASTA error in '%s': missing header line\n", filename);
        exit(EXIT_FAILURE);
    }

    if (length == 0) {
        free(bases);
        fprintf(stderr, "FASTA error in '%s': empty sequence\n", filename);
        exit(EXIT_FAILURE);
    }

    Strand *strand = strand_create(length);
    for (size_t i = 0; i < length; i++)
        strand->bases[i] = bases[i];

    free(bases);
    return strand;
}

