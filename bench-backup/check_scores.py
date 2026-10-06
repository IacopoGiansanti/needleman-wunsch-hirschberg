#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict

input_path = sys.argv[1] if len(sys.argv) > 1 else "results/results.csv"

instances = defaultdict(dict)
invalid = []
failed = []

with open(input_path, newline="") as f:
    for row in csv.DictReader(f):
        status = row["status"]
        if status != "ok":
            failed.append(row)
            continue

        if row["valid"] != "1":
            invalid.append(row)
            continue

        key = (
            int(row["n"]),
            int(row["m"]),
            int(row["seed"]),
            row["input_mode"],
            float(row["mutation_rate"]),
        )
        instances[key][row["algorithm"]] = int(row["score"])

mismatches = []
paired = 0

for key, scores in instances.items():
    if "nw" in scores and "hirschberg" in scores:
        paired += 1
        if scores["nw"] != scores["hirschberg"]:
            mismatches.append((key, scores))

if invalid:
    print(f"ERRORE: {len(invalid)} allineamenti non validi.")

if mismatches:
    print(f"ERRORE: {len(mismatches)} istanze con score ottimi diversi:")
    for key, scores in mismatches[:20]:
        print(" ", key, scores)

print(f"Istanze confrontabili NW/Hirschberg: {paired}")
print(f"Esecuzioni fallite/OOM registrate  : {len(failed)}")

if failed:
    by_alg = defaultdict(int)
    for row in failed:
        by_alg[row["algorithm"]] += 1
    print("Failure per algoritmo:", dict(by_alg))

if invalid or mismatches:
    sys.exit(1)

print("OK: tutti gli score confrontabili coincidono e gli allineamenti sono validi.")

