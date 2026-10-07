from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ------------------------------------------------------------
# Percorsi
# ------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "bench" / "results"

# Prende automaticamente l'ultima campagna "all"
campaigns = sorted(
    p for p in RESULTS_DIR.iterdir()
    if p.is_dir() and p.name.endswith("_all")
)

if not campaigns:
    raise RuntimeError("Nessuna campagna all trovata")

campaign = campaigns[-1]
summary_path = campaign / "summary.csv"

# Cartella in cui vengono salvati i grafici
output_dir = ROOT / "tesi-tex" / "figure"
output_dir.mkdir(exist_ok=True)


# ------------------------------------------------------------
# Caricamento dati
# ------------------------------------------------------------

df = pd.read_csv(summary_path)

# Utilizziamo solo:
# - sequenze sintetiche random
# - istanze quadrate n = m
random_square = df[
    (df["dataset_id"] == "synthetic_random")
    & (df["n"] == df["m"])
].copy()

random_square = random_square.sort_values("n")


# ------------------------------------------------------------
# Figura 1: tempo
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(6.4, 4.5))

for algorithm, label in [
    ("nw", "Needleman--Wunsch"),
    ("hirschberg", "Hirschberg"),
]:
    data = random_square[
        random_square["algorithm"] == algorithm
    ]

    ax.loglog(
        data["n"],
        data["wall_median_s"],
        marker="o",
        linewidth=1.5,
        label=label,
    )

ax.set_xlabel(r"Lunghezza delle sequenze $n=m$")
ax.set_ylabel("Wall time mediano [s]")

ax.grid(True, which="both", alpha=0.3)
ax.legend()

fig.tight_layout()

fig.savefig(
    output_dir / "scalabilita_tempo.pdf",
    bbox_inches="tight",
)

plt.close(fig)


# ------------------------------------------------------------
# Figura 2: memoria
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(6.4, 4.5))

for algorithm, label in [
    ("nw", "Needleman--Wunsch"),
    ("hirschberg", "Hirschberg"),
]:
    data = random_square[
        random_square["algorithm"] == algorithm
    ]

    ax.loglog(
        data["n"],
        data["rss_median_mib"],
        marker="o",
        linewidth=1.5,
        label=label,
    )

ax.set_xlabel(r"Lunghezza delle sequenze $n=m$")
ax.set_ylabel("Massimo RSS mediano [MiB]")

ax.grid(True, which="both", alpha=0.3)
ax.legend()

fig.tight_layout()

fig.savefig(
    output_dir / "scalabilita_memoria.pdf",
    bbox_inches="tight",
)

plt.close(fig)


# ------------------------------------------------------------
# Regressioni log-log utili per l'analisi
# ------------------------------------------------------------

print(f"Campagna: {campaign.name}")

for algorithm in ["nw", "hirschberg"]:
    data = random_square[
        random_square["algorithm"] == algorithm
    ]

    beta, alpha = np.polyfit(
        np.log10(data["n"]),
        np.log10(data["wall_median_s"]),
        1,
    )

    print(
        f"Tempo {algorithm}: "
        f"beta = {beta:.3f}"
    )

# Per la memoria NW escludiamo le dimensioni più piccole,
# dove l'overhead del processo pesa molto sul RSS.
nw_memory = random_square[
    (random_square["algorithm"] == "nw")
    & (random_square["n"] >= 1000)
]

beta, alpha = np.polyfit(
    np.log10(nw_memory["n"]),
    np.log10(nw_memory["rss_median_mib"]),
    1,
)

print(
    f"Memoria NW (n >= 1000): "
    f"beta = {beta:.3f}"
)
