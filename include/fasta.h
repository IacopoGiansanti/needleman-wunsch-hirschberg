#ifndef FASTA_H
#define FASTA_H

#include "strand.h"

/*
 * Legge un file FASTA contenente una singola sequenza di DNA e restituisce
 * uno Strand allocato dinamicamente.
 *
 * Il parser:
 *   - richiede una riga di intestazione che inizi con '>';
 *   - accetta sequenze distribuite su piu' righe;
 *   - ignora spazi bianchi nelle righe di sequenza;
 *   - accetta A, C, G, T sia maiuscole sia minuscole;
 *   - rifiuta caratteri non appartenenti all'alfabeto supportato;
 *   - rifiuta file FASTA con piu' di un record, per evitare ambiguita'.
 *
 * In caso di errore stampa un messaggio su stderr e termina il programma.
 * Il chiamante diventa proprietario dello Strand restituito e deve liberarlo
 * con strand_free().
 */
Strand *fasta_read_strand(const char *filename);

#endif // !FASTA_H
