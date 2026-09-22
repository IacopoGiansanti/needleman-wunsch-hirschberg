# Global Sequence Alignment: Needleman-Wunsch vs Hirschberg

Implementation in C of two algorithms for global sequence alignment:

- Needleman-Wunsch
- Hirschberg linear-space alignment

The project compares the two algorithms in terms of:

- alignment score
- execution time
- memory usage

## Scoring model

Default scoring:

- match: +1
- mismatch: -1
- gap: -1

## Build

```bash
make
=======
# needleman-wunsch-hirschberg
Implementation of Hirschberg and Needleman-Wunsch algorithm for the Global Alignment Problem
