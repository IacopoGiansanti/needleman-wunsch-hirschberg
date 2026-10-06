#!/usr/bin/env python3
"""Statistiche delle sole misure valide, con conteggio degli stati esclusi."""
import csv
import statistics
import sys
from collections import Counter, defaultdict
from check_scores import assess


def summarize(rows, output):
    # Individuare score discordanti anche se questo script è invocato da solo.
    bad_pairs = set()
    pairs = defaultdict(list)
    seen = set()
    for r in rows:
        key = r.get('pair_id') or (r.get('dataset_id', ''), r['n'], r['m'], r['seed'],
                r['input_mode'], r['mutation_rate'], r.get('rep', ''), r.get('orientation', ''))
        if (key, r['algorithm']) in seen:
            raise ValueError('CSV con duplicati: controllare prima di aggregare')
        seen.add((key, r['algorithm']))
        pairs[key].append(r)
    for key, pair in pairs.items():
        valid = [r for r in pair if r['status'] == 'ok' and r['valid'] == '1']
        if len(valid) == 2 and int(valid[0]['score']) != int(valid[1]['score']):
            bad_pairs.add(key)
    groups = defaultdict(list)
    for r in rows:
        key = (r.get('dataset_id', ''), r['algorithm'], r['n'], r['m'],
               r['input_mode'], r['mutation_rate'], r.get('orientation', ''))
        groups[key].append(r)
    keys = ['dataset_id', 'algorithm', 'n', 'm', 'input_mode', 'mutation_rate', 'orientation']
    fields = keys + ['attempts', 'runs', 'excluded', 'skipped_memory_budget',
        'score_pairs_matched', 'status_counts', 'wall_median_s', 'wall_mean_s', 'wall_stdev_s',
        'cpu_median_s', 'cpu_mean_s', 'cpu_stdev_s', 'rss_median_kb', 'rss_median_mib',
        'rss_mean_kb', 'rss_stdev_kb', 'score_min', 'score_max']
    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for key, group in sorted(groups.items()):
            valid = []
            for r in group:
                pk = r.get('pair_id') or (r.get('dataset_id', ''), r['n'], r['m'], r['seed'],
                     r['input_mode'], r['mutation_rate'], r.get('rep', ''), r.get('orientation', ''))
                if (r['status'] == 'ok' and r['valid'] == '1' and
                    r.get('score_check') != 'mismatch' and pk not in bad_pairs):
                    valid.append(r)
            counts = Counter(r['status'] for r in group)
            result = dict(zip(keys, key))
            result.update(attempts=len(group), runs=len(valid), excluded=len(group)-len(valid),
                skipped_memory_budget=counts['skipped_memory_budget'],
                score_pairs_matched=sum(r.get('score_check') == 'matched' for r in valid),
                status_counts=';'.join(f'{s}:{c}' for s, c in sorted(counts.items())))
            if valid:
                for column, prefix, unit in [('wall_time_s', 'wall', 's'),
                    ('cpu_time_s', 'cpu', 's'), ('max_rss_kb', 'rss', 'kb')]:
                    values = [float(r[column]) for r in valid]
                    result[f'{prefix}_median_{unit}'] = statistics.median(values)
                    result[f'{prefix}_mean_{unit}'] = statistics.mean(values)
                    result[f'{prefix}_stdev_{unit}'] = statistics.stdev(values) if len(values) > 1 else ''
                result['rss_median_mib'] = result['rss_median_kb'] / 1024
                scores = [int(r['score']) for r in valid]
                result.update(score_min=min(scores), score_max=max(scores))
            writer.writerow(result)
    print(f'Creato {output} ({len(groups)} gruppi). Deviazione standard campionaria; vuota con <2 misure.')


if __name__ == '__main__':
    source = sys.argv[1] if len(sys.argv) > 1 else 'results/results.csv'
    output = sys.argv[2] if len(sys.argv) > 2 else 'results/summary.csv'
    with open(source, newline='') as f:
        rows = list(csv.DictReader(f))
    # Rifiutare metadati corrotti; fallimenti registrati rimangono aggregabili.
    issues, _, _, _ = assess(rows)
    structural = [x for x in issues if 'metadati diversi' in x or 'duplicato' in x or 'sconosciuto' in x]
    if structural:
        sys.exit('Errore: ' + '; '.join(structural[:5]))
    summarize(rows, output)
