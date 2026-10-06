#!/usr/bin/env python3

import csv
import statistics
import sys
from collections import defaultdict

input_path = sys.argv[1] if len(sys.argv) > 1 else "results/results.csv"
output_path = sys.argv[2] if len(sys.argv) > 2 else "results/summary.csv"

groups = defaultdict(lambda: {"wall": [], "cpu": [], "rss": [], "scores": []})

with open(input_path, newline="") as f:
    for row in csv.DictReader(f):
        if row["status"] != "ok" or row["valid"] != "1":
            continue

        key = (
            row["algorithm"],
            int(row["n"]),
            int(row["m"]),
            row["input_mode"],
            float(row["mutation_rate"]),
        )

        d = groups[key]
        d["wall"].append(float(row["wall_time_s"]))
        d["cpu"].append(float(row["cpu_time_s"]))
        d["rss"].append(float(row["max_rss_kb"]))
        d["scores"].append(int(row["score"]))


def stdev(values):
    return statistics.stdev(values) if len(values) > 1 else 0.0


fieldnames = [
    "algorithm", "n", "m", "input_mode", "mutation_rate", "runs",
    "wall_median_s", "wall_mean_s", "wall_stdev_s",
    "cpu_median_s", "cpu_mean_s", "cpu_stdev_s",
    "rss_median_kb", "rss_median_mib", "rss_mean_kb", "rss_stdev_kb",
    "score_min", "score_max",
]

with open(output_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    for key in sorted(groups, key=lambda x: (x[3], x[1], x[2], x[0])):
        algorithm, n, m, mode, mutation_rate = key
        d = groups[key]
        rss_median = statistics.median(d["rss"])

        writer.writerow({
            "algorithm": algorithm,
            "n": n,
            "m": m,
            "input_mode": mode,
            "mutation_rate": mutation_rate,
            "runs": len(d["wall"]),
            "wall_median_s": statistics.median(d["wall"]),
            "wall_mean_s": statistics.mean(d["wall"]),
            "wall_stdev_s": stdev(d["wall"]),
            "cpu_median_s": statistics.median(d["cpu"]),
            "cpu_mean_s": statistics.mean(d["cpu"]),
            "cpu_stdev_s": stdev(d["cpu"]),
            "rss_median_kb": rss_median,
            "rss_median_mib": rss_median / 1024.0,
            "rss_mean_kb": statistics.mean(d["rss"]),
            "rss_stdev_kb": stdev(d["rss"]),
            "score_min": min(d["scores"]),
            "score_max": max(d["scores"]),
        })

print(f"Creato {output_path} con {len(groups)} gruppi sperimentali.")

