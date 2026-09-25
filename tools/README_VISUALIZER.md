# Alignment visualizer

This directory contains a visualization layer for the project without modifying the
alignment algorithms in `src/`.

The visualizer uses the real C implementations as the source of truth for the final
alignment. The small helper `alignment_dump.c` links against the existing project
sources and exposes the returned `Alignment` in a machine-readable form. Python then
reconstructs only the information needed for figures and animations.

## What it shows

### Needleman-Wunsch

- score-matrix construction;
- explored states in the edit/alignment graph;
- traceback from sink to source;
- final optimal path;
- match/mismatch/insertion/deletion labels on the path;
- final alignment returned by the real C implementation.

### Hirschberg

- active divide-and-conquer rectangle at every recursive call;
- middle vertex and middle edge selected by the same rules used in `src/Hirschberg.c`
  and `src/middle_edge.c`;
- prefix/suffix score vectors for the current split;
- accumulation of middle edges;
- final path reconstructed from the alignment returned by the real C implementation.

The complete graph is drawn only for explanatory purposes. It is **not** used by the
actual Hirschberg implementation and must not be included in memory benchmarks.

## Dependencies

The C helper only needs the same compiler already used by the project.

Python dependencies:

```bash
python3 -m pip install -r tools/requirements.txt
```

## Build the C bridge

The visualizer builds it automatically when needed, but it can also be compiled with:

```bash
make -C tools
```

This creates:

```text
tools/alignment_dump
```

No source file under `src/` or `include/` is changed.

## Basic usage

From the repository root:

```bash
python3 tools/visualizer.py fasta/HIV-1-95TNIH047.fasta \
                            fasta/HIV-1-95TNIH022.fasta \
                            --algorithm both \
                            --slice1 0:12 \
                            --slice2 0:12 \
                            --name hiv_demo
```

Generated files are written by default to `visualization/`:

```text
visualization/hiv_demo_nw_final.png
visualization/hiv_demo_nw_final.svg
visualization/hiv_demo_nw.gif
visualization/hiv_demo_hirschberg_final.png
visualization/hiv_demo_hirschberg_final.svg
visualization/hiv_demo_hirschberg.gif
```

## Why slices are recommended

The real HIV/SARS-CoV-2 FASTA files contain thousands of nucleotides. Their edit graph
would contain millions of nodes/edges and would not be readable as a figure.

Use half-open, 0-based slices:

```bash
--slice1 100:120 --slice2 105:125
```

The actual C algorithms are run on exactly those subsequences. This means the plotted
optimal path is still the output of the project implementation, just on a deliberately
small instance suitable for visualization.

The visualizer refuses grids with a side longer than 35 by default. `--force` overrides
this limit, but very large plots are not useful in practice.

## Scoring

Defaults match the current project:

```text
match     +1
mismatch  -1
gap       -1
```

They can be changed for a visualization:

```bash
python3 tools/visualizer.py a.fasta b.fasta \
    --algorithm nw \
    --match 2 --mismatch -1 --gap -2
```

The same values are passed to the real C implementation and to the visualization
simulation.

## Static figures only

To skip GIF creation:

```bash
python3 tools/visualizer.py a.fasta b.fasta --algorithm both --no-animation
```

## Animation controls

```bash
--fps 4
--max-frames 160
```

For long-but-still-readable grids, `--max-frames` samples the DP construction so the GIF
does not contain hundreds or thousands of almost identical frames. The final C result is
unchanged.

## Output interpretation

Path labels are:

```text
M  match
X  mismatch
I  insertion (gap in sequence 1)
D  deletion  (gap in sequence 2)
```

For Needleman-Wunsch, the DP simulator reproduces the recurrence and traceback
priority in `src/Needleman-Wunsch.c` (`DIAGONAL > UP > LEFT` when ties are present).
The final path is nevertheless always reconstructed from the alignment returned by the
compiled C implementation.

For Hirschberg, the visualizer mirrors the project's current middle-node and middle-edge
tie rules. The orange dashed edges are recursive middle edges; the red path is the final
alignment returned by the C implementation.

## Suggested use in the thesis/project report

A compact figure pair works well:

1. Needleman-Wunsch: complete DP matrix + optimal path.
2. Hirschberg: same conceptual edit graph with the recursive middle edges highlighted.

This visually communicates the key distinction: both algorithms return an optimal global
alignment, while Hirschberg avoids retaining the complete quadratic DP/backtracking
structure during the actual computation.
