#!/usr/bin/env python3
"""Visualize the alignment graph used by this repository.

This tool is intentionally external to the implementation under src/.
It calls the real C implementations through tools/alignment_dump and uses
Python only to reconstruct visualization data and animations.

Outputs:
  * final PNG
  * final SVG
  * animated GIF (optional)

For Needleman-Wunsch the animation shows matrix construction followed by
traceback. For Hirschberg the animation shows the recursive subproblems and
middle edges, then reveals the final path returned by the real C code.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

Point = Tuple[int, int]


@dataclass(frozen=True)
class Scoring:
    match: int = 1
    mismatch: int = -1
    gap: int = -1

    def score(self, a: str, b: str) -> int:
        if a == "-" or b == "-":
            return self.gap
        return self.match if a == b else self.mismatch


@dataclass
class CAlignment:
    score: int
    seq1: str
    seq2: str

    @property
    def path(self) -> List[Point]:
        i = j = 0
        result: List[Point] = [(0, 0)]
        for a, b in zip(self.seq1, self.seq2):
            if a != "-":
                i += 1
            if b != "-":
                j += 1
            result.append((i, j))
        return result

    @property
    def operations(self) -> List[str]:
        ops: List[str] = []
        for a, b in zip(self.seq1, self.seq2):
            if a == "-":
                ops.append("I")
            elif b == "-":
                ops.append("D")
            elif a == b:
                ops.append("M")
            else:
                ops.append("X")
        return ops


@dataclass
class NWTrace:
    matrix: np.ndarray
    pred: List[List[str]]
    fill_order: List[Point]
    simulated_path: List[Point]


@dataclass(frozen=True)
class HirschEvent:
    depth: int
    active_rect: Tuple[int, int, int, int]  # top,bottom,left,right in ORIGINAL coordinates
    middle_vertex: Point
    edge_start: Point
    edge_end: Point
    prefix_scores: Tuple[int, ...]
    suffix_scores: Tuple[int, ...]


# ---------------------------------------------------------------------------
# Input / real C implementation bridge
# ---------------------------------------------------------------------------


def read_fasta(path: Path) -> str:
    header_seen = False
    sequence: List[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header_seen and sequence:
                    raise ValueError(f"{path}: multiple FASTA records are not supported")
                header_seen = True
                continue
            if not header_seen:
                raise ValueError(f"{path}: expected FASTA header beginning with '>'")
            for ch in line.upper():
                if ch not in "ACGT":
                    raise ValueError(f"{path}: unsupported DNA base {ch!r}")
                sequence.append(ch)
    if not header_seen:
        raise ValueError(f"{path}: missing FASTA header")
    if not sequence:
        raise ValueError(f"{path}: empty sequence")
    return "".join(sequence)


def parse_slice(spec: Optional[str], length: int) -> Tuple[int, int]:
    if spec is None:
        return (0, length)
    try:
        left_s, right_s = spec.split(":", 1)
        left = int(left_s) if left_s else 0
        right = int(right_s) if right_s else length
    except Exception as exc:
        raise ValueError(f"invalid slice {spec!r}; expected START:END") from exc
    if left < 0 or right < 0 or left > right or right > length:
        raise ValueError(f"slice {spec!r} is outside sequence length {length}")
    if left == right:
        raise ValueError("empty slices are not useful for this visualizer")
    return left, right


def ensure_dump_binary(repo_root: Path) -> Path:
    tool_dir = repo_root / "tools"
    binary = tool_dir / "alignment_dump"
    sources = [
        tool_dir / "alignment_dump.c",
        repo_root / "src" / "Needleman-Wunsch.c",
        repo_root / "src" / "Hirschberg.c",
        repo_root / "src" / "middle_edge.c",
    ]
    needs_build = not binary.exists()
    if binary.exists():
        binary_mtime = binary.stat().st_mtime
        needs_build = any(p.exists() and p.stat().st_mtime > binary_mtime for p in sources)
    if needs_build:
        subprocess.run(["make", "-C", str(tool_dir)], check=True)
    return binary


def write_temp_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w", encoding="ascii") as fh:
        fh.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            fh.write(sequence[start:start + 80] + "\n")


def run_real_alignment(repo_root: Path,
                       algorithm: str,
                       seq1: str,
                       seq2: str,
                       scoring: Scoring) -> CAlignment:
    binary = ensure_dump_binary(repo_root)
    with tempfile.TemporaryDirectory(prefix="alignment-vis-") as tmp:
        tmpdir = Path(tmp)
        f1 = tmpdir / "seq1.fasta"
        f2 = tmpdir / "seq2.fasta"
        write_temp_fasta(f1, "visualizer_seq1", seq1)
        write_temp_fasta(f2, "visualizer_seq2", seq2)
        proc = subprocess.run(
            [str(binary), algorithm, str(f1), str(f2),
             str(scoring.match), str(scoring.mismatch), str(scoring.gap)],
            check=True,
            text=True,
            capture_output=True,
        )

    values = {}
    for line in proc.stdout.splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            values[key] = value.strip()
    try:
        result = CAlignment(
            score=int(values["SCORE"]),
            seq1=values["SEQ1"],
            seq2=values["SEQ2"],
        )
    except KeyError as exc:
        raise RuntimeError(f"unexpected alignment_dump output:\n{proc.stdout}") from exc
    if len(result.seq1) != len(result.seq2):
        raise RuntimeError("real C implementation returned rows of different length")
    return result


# ---------------------------------------------------------------------------
# Faithful visualization simulations
# ---------------------------------------------------------------------------


def nw_trace(seq1: str, seq2: str, scoring: Scoring) -> NWTrace:
    """Replicate the recurrence/tie-breaking of src/Needleman-Wunsch.c."""
    n, m = len(seq1), len(seq2)
    score = np.zeros((n + 1, m + 1), dtype=int)
    pred = [["" for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + scoring.gap
        pred[i][0] = "U"
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + scoring.gap
        pred[0][j] = "L"

    fill_order: List[Point] = []
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            up = int(score[i - 1, j]) + scoring.gap
            left = int(score[i, j - 1]) + scoring.gap
            diag = int(score[i - 1, j - 1]) + scoring.score(seq1[i - 1], seq2[j - 1])
            best = max(up, left, diag)
            score[i, j] = best
            moves = ""
            if best == up:
                moves += "U"
            if best == left:
                moves += "L"
            if best == diag:
                moves += "D"
            pred[i][j] = moves
            fill_order.append((i, j))

    # C traceback priority is DIAGONAL > UP > LEFT.
    i, j = n, m
    reverse_path = [(i, j)]
    while i > 0 or j > 0:
        moves = pred[i][j]
        if "D" in moves:
            i -= 1
            j -= 1
        elif "U" in moves:
            i -= 1
        else:
            j -= 1
        reverse_path.append((i, j))
    simulated_path = list(reversed(reverse_path))
    return NWTrace(score, pred, fill_order, simulated_path)


def _nw_score_columns(rows: str, r0: int, r1: int,
                      cols: str, c0: int, c1: int,
                      scoring: Scoring) -> List[int]:
    n = r1 - r0
    m = c1 - c0
    previous = [0] * (n + 1)
    current = [0] * (n + 1)
    for i in range(1, n + 1):
        previous[i] = previous[i - 1] + scoring.gap
    for j in range(1, m + 1):
        current[0] = previous[0] + scoring.gap
        for i in range(1, n + 1):
            up = current[i - 1] + scoring.gap
            left = previous[i] + scoring.gap
            diagonal = previous[i - 1] + scoring.score(rows[r0 + i - 1], cols[c0 + j - 1])
            current[i] = max(up, left, diagonal)
        previous, current = current, previous
    return previous


def _nw_score_reverse(rows: str, r0: int, r1: int,
                      cols: str, c0: int, c1: int,
                      scoring: Scoring) -> List[int]:
    n = r1 - r0
    m = c1 - c0
    previous = [0] * (n + 1)
    current = [0] * (n + 1)
    previous[n] = 0
    for i in range(n - 1, -1, -1):
        previous[i] = previous[i + 1] + scoring.gap
    for j in range(m - 1, -1, -1):
        current[n] = previous[n] + scoring.gap
        for i in range(n - 1, -1, -1):
            down = current[i + 1] + scoring.gap
            right = previous[i] + scoring.gap
            diagonal = previous[i + 1] + scoring.score(rows[r0 + i], cols[c0 + j])
            current[i] = max(down, right, diagonal)
        previous, current = current, previous
    return previous


def hirschberg_events(seq1: str, seq2: str, scoring: Scoring) -> List[HirschEvent]:
    """Mirror the middle-edge choices in src/Hirschberg.c for visualization.

    The original implementation makes the shorter sequence the row dimension.
    Events are mapped back to the original (seq1, seq2) edit-graph coordinates.
    """
    swapped = len(seq1) > len(seq2)
    rows = seq2 if swapped else seq1
    cols = seq1 if swapped else seq2
    events: List[HirschEvent] = []

    def map_point(p: Point) -> Point:
        r, c = p
        return (c, r) if swapped else (r, c)

    def map_rect(r0: int, r1: int, c0: int, c1: int) -> Tuple[int, int, int, int]:
        if swapped:
            return (c0, c1, r0, r1)
        return (r0, r1, c0, c1)

    def rec(r0: int, r1: int, c0: int, c1: int, depth: int) -> None:
        if c0 == c1 or r0 == r1:
            return

        n = r1 - r0
        m = c1 - c0
        mid_j = c0 + m // 2

        prefix = _nw_score_columns(rows, r0, r1, cols, c0, mid_j, scoring)
        suffix = _nw_score_reverse(rows, r0, r1, cols, mid_j, c1, scoring)

        mid_offset = 0
        best = prefix[0] + suffix[0]
        for i in range(1, n + 1):
            candidate = prefix[i] + suffix[i]
            if best < candidate:  # strict, matching C: first maximum wins
                best = candidate
                mid_offset = i

        mid_i = r0 + mid_offset
        suffix_next = _nw_score_reverse(rows, r0, r1, cols, mid_j + 1, c1, scoring)

        local_i = mid_i - r0
        edge_start = (mid_i, mid_j)
        if mid_i == r1:
            edge_end = (mid_i, mid_j + 1)
        else:
            down = suffix[local_i + 1] + scoring.gap
            right = suffix_next[local_i] + scoring.gap
            diagonal = suffix_next[local_i + 1] + scoring.score(rows[mid_i], cols[mid_j])
            final = max(down, right, diagonal)
            # mid_edge.c priority: diagonal > right > down
            if final == diagonal:
                edge_end = (mid_i + 1, mid_j + 1)
            elif final == right:
                edge_end = (mid_i, mid_j + 1)
            else:
                edge_end = (mid_i + 1, mid_j)

        events.append(HirschEvent(
            depth=depth,
            active_rect=map_rect(r0, r1, c0, c1),
            middle_vertex=map_point((mid_i, mid_j)),
            edge_start=map_point(edge_start),
            edge_end=map_point(edge_end),
            prefix_scores=tuple(prefix),
            suffix_scores=tuple(suffix),
        ))

        rec(r0, mid_i, c0, mid_j, depth + 1)
        rec(edge_end[0], r1, edge_end[1], c1, depth + 1)

    rec(0, len(rows), 0, len(cols), 0)
    return events


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------


def operation_for_edge(a: Point, b: Point, seq1: str, seq2: str) -> str:
    di = b[0] - a[0]
    dj = b[1] - a[1]
    if di == 1 and dj == 1:
        return "M" if seq1[a[0]] == seq2[a[1]] else "X"
    if di == 1:
        return "D"
    return "I"


def draw_graph(ax,
               seq1: str,
               seq2: str,
               path: Optional[Sequence[Point]] = None,
               path_limit: Optional[int] = None,
               active_rect: Optional[Tuple[int, int, int, int]] = None,
               middle_edges: Optional[Sequence[Tuple[Point, Point]]] = None,
               current_edge: Optional[Tuple[Point, Point]] = None,
               current_vertex: Optional[Point] = None,
               explored: Optional[Iterable[Point]] = None,
               title: str = "Alignment graph") -> None:
    n, m = len(seq1), len(seq2)

    def xy(p: Point) -> Tuple[float, float]:
        i, j = p
        return float(j), float(n - i)

    # Full conceptual graph in a subtle background.
    for i in range(n + 1):
        for j in range(m + 1):
            x, y = xy((i, j))
            if i < n:
                x2, y2 = xy((i + 1, j))
                ax.plot([x, x2], [y, y2], color="#c8c8c8", linewidth=0.55, alpha=0.45, zorder=1)
            if j < m:
                x2, y2 = xy((i, j + 1))
                ax.plot([x, x2], [y, y2], color="#c8c8c8", linewidth=0.55, alpha=0.45, zorder=1)
            if i < n and j < m:
                x2, y2 = xy((i + 1, j + 1))
                ax.plot([x, x2], [y, y2], color="#d8d8d8", linewidth=0.55, alpha=0.40, zorder=1)

    if explored is not None:
        pts = list(explored)
        if pts:
            ax.scatter([xy(p)[0] for p in pts], [xy(p)[1] for p in pts],
                       s=21, color="#6baed6", alpha=0.75, zorder=2)

    # All nodes.
    all_pts = [(i, j) for i in range(n + 1) for j in range(m + 1)]
    ax.scatter([xy(p)[0] for p in all_pts], [xy(p)[1] for p in all_pts],
               s=13, color="#8c8c8c", alpha=0.45, zorder=2)

    if active_rect is not None:
        top, bottom, left, right = active_rect
        x0, y_top = xy((top, left))
        x1, y_bottom = xy((bottom, right))
        ax.add_patch(Rectangle((x0, y_bottom), x1 - x0, y_top - y_bottom,
                               fill=False, edgecolor="#ff8c00", linewidth=2.4, zorder=4))

    if middle_edges:
        for a, b in middle_edges:
            x1, y1 = xy(a)
            x2, y2 = xy(b)
            ax.plot([x1, x2], [y1, y2], color="#ff8c00", linestyle="--", linewidth=2.2, zorder=5)

    if current_edge:
        a, b = current_edge
        x1, y1 = xy(a)
        x2, y2 = xy(b)
        ax.plot([x1, x2], [y1, y2], color="#ffb000", linewidth=4.0, zorder=7)

    if path:
        shown = list(path if path_limit is None else path[:path_limit])
        if len(shown) >= 2:
            for a, b in zip(shown, shown[1:]):
                x1, y1 = xy(a)
                x2, y2 = xy(b)
                ax.plot([x1, x2], [y1, y2], color="#c62828", linewidth=3.0, zorder=8)
                xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
                ax.text(xm, ym, operation_for_edge(a, b, seq1, seq2),
                        fontsize=7, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="#444444", alpha=0.88), zorder=9)
        if shown:
            ax.scatter([xy(p)[0] for p in shown], [xy(p)[1] for p in shown],
                       s=28, color="#c62828", zorder=9)

    if current_vertex:
        x, y = xy(current_vertex)
        ax.scatter([x], [y], s=80, color="#ffb000", edgecolor="black", linewidth=1.0, zorder=10)

    for j, ch in enumerate(seq2, 1):
        ax.text(j, n + 0.55, ch, ha="center", va="center", fontweight="bold")
    ax.text(0, n + 0.55, "-", ha="center", va="center", fontweight="bold")
    for i, ch in enumerate(seq1, 1):
        ax.text(-0.55, n - i, ch, ha="center", va="center", fontweight="bold")
    ax.text(-0.55, n, "-", ha="center", va="center", fontweight="bold")

    ax.set_xlim(-1.0, m + 0.7)
    ax.set_ylim(-0.7, n + 1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title)


def draw_matrix(ax,
                seq1: str,
                seq2: str,
                trace: NWTrace,
                visible_fill: int,
                path_cells: Optional[Sequence[Point]] = None,
                current: Optional[Point] = None,
                title: str = "Needleman-Wunsch score matrix") -> None:
    n, m = len(seq1), len(seq2)
    shown = np.full(trace.matrix.shape, np.nan, dtype=float)
    shown[:, 0] = trace.matrix[:, 0]
    shown[0, :] = trace.matrix[0, :]
    for p in trace.fill_order[:visible_fill]:
        shown[p] = trace.matrix[p]

    ax.imshow(np.ma.masked_invalid(shown), origin="upper", cmap="viridis", alpha=0.55)
    ax.set_xticks(range(m + 1), ["-"] + list(seq2))
    ax.set_yticks(range(n + 1), ["-"] + list(seq1))
    ax.set_xticks(np.arange(-.5, m + 1, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n + 1, 1), minor=True)
    ax.grid(which="minor", linewidth=0.6, alpha=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(n + 1):
        for j in range(m + 1):
            if not math.isnan(shown[i, j]):
                ax.text(j, i, str(int(shown[i, j])), ha="center", va="center", fontsize=8)

    if path_cells:
        for i, j in path_cells:
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1,
                                   fill=False, edgecolor="#c62828", linewidth=2.0, zorder=5))
    if current:
        i, j = current
        ax.add_patch(Rectangle((j - .5, i - .5), 1, 1,
                               fill=False, edgecolor="#ffb000", linewidth=3.2, zorder=6))
    ax.set_title(title)


def add_alignment_text(fig, result: CAlignment, algorithm_label: str) -> None:
    middle = []
    for a, b in zip(result.seq1, result.seq2):
        if a == "-" or b == "-":
            middle.append(" ")
        elif a == b:
            middle.append("|")
        else:
            middle.append(".")
    fig.text(
        0.5, 0.015,
        f"{algorithm_label}  score={result.score}\n"
        f"seq1  {result.seq1}\n"
        f"      {''.join(middle)}\n"
        f"seq2  {result.seq2}",
        ha="center", va="bottom", family="monospace", fontsize=9,
    )


def save_nw_static(seq1: str, seq2: str, trace: NWTrace, result: CAlignment,
                   png: Path, svg: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    draw_matrix(axes[0], seq1, seq2, trace, len(trace.fill_order), result.path,
                title="Needleman-Wunsch: final DP matrix")
    draw_graph(axes[1], seq1, seq2, path=result.path,
               title="Edit graph and optimal path")
    add_alignment_text(fig, result, "Needleman-Wunsch (real C implementation)")
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def save_hirschberg_static(seq1: str, seq2: str, events: List[HirschEvent], result: CAlignment,
                           png: Path, svg: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    middle_edges = [(e.edge_start, e.edge_end) for e in events]
    draw_graph(axes[0], seq1, seq2, path=result.path, middle_edges=middle_edges,
               title="Hirschberg: middle edges + final optimal path")
    axes[1].axis("off")
    axes[1].set_title("Divide-and-conquer events")
    if events:
        lines = ["#  depth   rectangle [rows x cols]    middle edge"]
        for k, e in enumerate(events[:28], 1):
            t, b, l, r = e.active_rect
            lines.append(
                f"{k:>2}  {e.depth:>5}   [{t}:{b}] x [{l}:{r}]   "
                f"{e.edge_start}->{e.edge_end}"
            )
        if len(events) > 28:
            lines.append(f"... {len(events) - 28} more events ...")
        axes[1].text(0.02, 0.97, "\n".join(lines), va="top", ha="left",
                     family="monospace", fontsize=8.5, transform=axes[1].transAxes)
    else:
        axes[1].text(0.5, 0.5, "No recursive split was necessary.", ha="center")
    add_alignment_text(fig, result, "Hirschberg (real C implementation)")
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def sampled_counts(total: int, max_frames: int) -> List[int]:
    if total <= max_frames:
        return list(range(1, total + 1))
    raw = np.linspace(1, total, max_frames, dtype=int)
    return sorted(set(int(x) for x in raw))


def animate_nw(seq1: str, seq2: str, trace: NWTrace, result: CAlignment,
               out: Path, fps: int, max_frames: int) -> None:
    fill_counts = sampled_counts(len(trace.fill_order), max(1, max_frames - len(result.path)))
    # traceback is shown from sink back to source, like the C implementation.
    traceback_path = list(reversed(result.path))
    trace_counts = sampled_counts(len(traceback_path), min(len(traceback_path), max_frames // 2 + 1))
    frames = [("fill", count) for count in fill_counts] + [("trace", count) for count in trace_counts]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    def update(frame):
        for ax in axes:
            ax.clear()
        phase, count = frame
        if phase == "fill":
            current = trace.fill_order[count - 1] if trace.fill_order else (0, 0)
            explored = {(i, 0) for i in range(len(seq1) + 1)} | {(0, j) for j in range(len(seq2) + 1)}
            explored |= set(trace.fill_order[:count])
            draw_matrix(axes[0], seq1, seq2, trace, count, current=current,
                        title=f"DP construction: {count}/{len(trace.fill_order)} cells")
            draw_graph(axes[1], seq1, seq2, explored=explored, current_vertex=current,
                       title="Edit graph: explored DP states")
        else:
            shown_back = traceback_path[:count]
            current = shown_back[-1]
            draw_matrix(axes[0], seq1, seq2, trace, len(trace.fill_order),
                        path_cells=shown_back, current=current,
                        title=f"Traceback: {count}/{len(traceback_path)} vertices")
            # reverse again so path edges display in source->sink direction.
            shown_forward = list(reversed(shown_back))
            draw_graph(axes[1], seq1, seq2, path=shown_forward, current_vertex=current,
                       title="Optimal path being recovered")
    add_alignment_text(fig, result, "Needleman-Wunsch (real C result)")
    fig.tight_layout(rect=(0, 0.13, 1, 1))

    anim = FuncAnimation(fig, update, frames=frames, interval=1000 / max(1, fps), repeat=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)


def animate_hirschberg(seq1: str, seq2: str, events: List[HirschEvent], result: CAlignment,
                       out: Path, fps: int, max_frames: int) -> None:
    path_counts = sampled_counts(len(result.path), min(max_frames // 2 + 1, len(result.path)))
    event_indices = list(range(len(events)))
    if len(event_indices) > max_frames // 2:
        event_indices = [x - 1 for x in sampled_counts(len(events), max_frames // 2)]
    frames = [("event", k) for k in event_indices] + [("path", c) for c in path_counts]
    if not frames:
        frames = [("path", len(result.path))]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    def update(frame):
        for ax in axes:
            ax.clear()
        phase, value = frame
        if phase == "event":
            idx = value
            e = events[idx]
            past = [(x.edge_start, x.edge_end) for x in events[:idx]]
            draw_graph(axes[0], seq1, seq2,
                       active_rect=e.active_rect,
                       middle_edges=past,
                       current_edge=(e.edge_start, e.edge_end),
                       current_vertex=e.middle_vertex,
                       title=f"Recursive split {idx + 1}/{len(events)} (depth {e.depth})")
            axes[1].axis("off")
            axes[1].set_title("Current Hirschberg subproblem")
            t, b, l, r = e.active_rect
            vector_len = min(len(e.prefix_scores), 18)
            prefix = list(e.prefix_scores[:vector_len])
            suffix = list(e.suffix_scores[:vector_len])
            detail = (
                f"depth          : {e.depth}\n"
                f"active rows    : [{t}, {b}]\n"
                f"active columns : [{l}, {r}]\n"
                f"middle vertex  : {e.middle_vertex}\n"
                f"middle edge    : {e.edge_start} -> {e.edge_end}\n\n"
                f"prefix scores  : {prefix}{' ...' if len(e.prefix_scores) > vector_len else ''}\n"
                f"suffix scores  : {suffix}{' ...' if len(e.suffix_scores) > vector_len else ''}\n\n"
                "Only linear-size score vectors are conceptually needed;\n"
                "the full graph shown at left exists only for visualization."
            )
            axes[1].text(0.02, 0.94, detail, va="top", ha="left",
                         family="monospace", fontsize=9, transform=axes[1].transAxes)
        else:
            count = value
            middle_edges = [(e.edge_start, e.edge_end) for e in events]
            draw_graph(axes[0], seq1, seq2, path=result.path, path_limit=count,
                       middle_edges=middle_edges,
                       current_vertex=result.path[min(count, len(result.path)) - 1],
                       title=f"Final path assembly: {count}/{len(result.path)} vertices")
            axes[1].axis("off")
            axes[1].set_title("Alignment returned by the C implementation")
            upto = max(0, min(count - 1, len(result.seq1)))
            partial1 = result.seq1[:upto]
            partial2 = result.seq2[:upto]
            ops = "".join(result.operations[:upto])
            axes[1].text(0.02, 0.75,
                         f"seq1: {partial1}\n"
                         f"ops : {ops}\n"
                         f"seq2: {partial2}\n\n"
                         "M=match, X=mismatch, I=insertion, D=deletion",
                         va="top", ha="left", family="monospace", fontsize=10,
                         transform=axes[1].transAxes)
    add_alignment_text(fig, result, "Hirschberg (real C result)")
    fig.tight_layout(rect=(0, 0.13, 1, 1))

    anim = FuncAnimation(fig, update, frames=frames, interval=1000 / max(1, fps), repeat=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Visualize the edit graph and algorithm evolution without modifying src/."
    )
    p.add_argument("fasta1", type=Path)
    p.add_argument("fasta2", type=Path)
    p.add_argument("--algorithm", "-a", choices=["nw", "hirschberg", "both"], default="both")
    p.add_argument("--slice1", metavar="START:END",
                   help="visualize only this 0-based half-open slice of FASTA 1")
    p.add_argument("--slice2", metavar="START:END",
                   help="visualize only this 0-based half-open slice of FASTA 2")
    p.add_argument("--match", type=int, default=1)
    p.add_argument("--mismatch", type=int, default=-1)
    p.add_argument("--gap", type=int, default=-1)
    p.add_argument("--output-dir", "-o", type=Path, default=Path("visualization"))
    p.add_argument("--name", default="alignment", help="output filename prefix")
    p.add_argument("--fps", type=int, default=3)
    p.add_argument("--max-frames", type=int, default=180,
                   help="cap GIF frames by sampling long animations")
    p.add_argument("--no-animation", action="store_true")
    p.add_argument("--max-side", type=int, default=35,
                   help="refuse larger visualization grids unless --force is used")
    p.add_argument("--force", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    repo_root = Path(__file__).resolve().parent.parent

    try:
        raw1 = read_fasta(args.fasta1)
        raw2 = read_fasta(args.fasta2)
        s1 = parse_slice(args.slice1, len(raw1))
        s2 = parse_slice(args.slice2, len(raw2))
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    seq1 = raw1[s1[0]:s1[1]]
    seq2 = raw2[s2[0]:s2[1]]

    if not args.force and (len(seq1) > args.max_side or len(seq2) > args.max_side):
        print(
            f"error: requested grid is {len(seq1)} x {len(seq2)}. "
            f"For a readable visualization keep each side <= {args.max_side}.\n"
            "Use --slice1/--slice2 (recommended for real FASTA files) or --force.",
            file=sys.stderr,
        )
        return 2

    scoring = Scoring(args.match, args.mismatch, args.gap)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    algorithms = ["nw", "hirschberg"] if args.algorithm == "both" else [args.algorithm]

    print("VISUALIZATION INPUT")
    print(f"  FASTA 1 : {args.fasta1} [{s1[0]}:{s1[1]}] -> {len(seq1)} nt")
    print(f"  FASTA 2 : {args.fasta2} [{s2[0]}:{s2[1]}] -> {len(seq2)} nt")
    print(f"  scoring : match={scoring.match:+d} mismatch={scoring.mismatch:+d} gap={scoring.gap:+d}")

    nw_sim = nw_trace(seq1, seq2, scoring)

    for algorithm in algorithms:
        result = run_real_alignment(repo_root, algorithm, seq1, seq2, scoring)
        print(f"\n{algorithm.upper()}")
        print(f"  real C score : {result.score}")
        print(f"  seq1         : {result.seq1}")
        print(f"  seq2         : {result.seq2}")

        prefix = args.output_dir / f"{args.name}_{algorithm}"
        png = prefix.with_name(prefix.name + "_final.png")
        svg = prefix.with_name(prefix.name + "_final.svg")
        gif = prefix.with_name(prefix.name + ".gif")

        if algorithm == "nw":
            if result.score != int(nw_sim.matrix[-1, -1]):
                raise RuntimeError("visual simulation score differs from C Needleman-Wunsch")
            if result.path != nw_sim.simulated_path:
                print("  note: C path differs from simulated NW tie-breaking; final path uses C output")
            save_nw_static(seq1, seq2, nw_sim, result, png, svg)
            if not args.no_animation:
                animate_nw(seq1, seq2, nw_sim, result, gif, args.fps, args.max_frames)
        else:
            events = hirschberg_events(seq1, seq2, scoring)
            save_hirschberg_static(seq1, seq2, events, result, png, svg)
            if not args.no_animation:
                animate_hirschberg(seq1, seq2, events, result, gif, args.fps, args.max_frames)
            print(f"  recursive middle-edge events : {len(events)}")

        print(f"  PNG : {png}")
        print(f"  SVG : {svg}")
        if not args.no_animation:
            print(f"  GIF : {gif}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
