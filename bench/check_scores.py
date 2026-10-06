#!/usr/bin/env python3
"""Verifica ogni ripetizione, senza sovrascrivere istanze o nascondere invalidi."""
import csv
import sys
from collections import Counter, defaultdict


def assess(rows):
    issues = []
    pairs = defaultdict(dict)
    statuses = Counter()
    for line, row in enumerate(rows, 2):
        status = row['status']
        statuses[status] += 1
        key = (row.get('pair_id') or (row.get('dataset_id', ''), row['n'], row['m'],
               row['seed'], row['input_mode'], row['mutation_rate'],
               row.get('rep', ''), row.get('orientation', '')))
        alg = row['algorithm']
        if alg not in ('nw', 'hirschberg'):
            issues.append(f'riga {line}: algoritmo sconosciuto')
        if alg in pairs[key]:
            issues.append(f'riga {line}: duplicato {key} / {alg}')
        pairs[key][alg] = row
        if status != 'ok' and status != 'skipped_memory_budget':
            issues.append(f'riga {line}: {status}')
        if status == 'ok' and row['valid'] != '1':
            issues.append(f'riga {line}: allineamento non valido')
        if status == 'score_mismatch' or row.get('score_check') == 'mismatch':
            issues.append(f'riga {line}: score discordanti')
    paired = unpaired = 0
    for key, pair in pairs.items():
        if set(pair) != {'nw', 'hirschberg'}:
            issues.append(f'{key}: risultato di un algoritmo mancante')
        good = [r for r in pair.values() if r['status'] == 'ok' and r['valid'] == '1']
        if len(pair) == 2:
            a, b = pair.values()
            for field in ('n', 'm', 'seed', 'input_mode', 'mutation_rate', 'dataset_id',
                          'rep', 'orientation', 'input1_sha256', 'input2_sha256'):
                if a.get(field, '') != b.get(field, ''):
                    issues.append(f'{key}: metadati diversi ({field})')
        if len(good) == 2:
            paired += 1
            if int(good[0]['score']) != int(good[1]['score']):
                issues.append(f'{key}: score diversi')
        else:
            unpaired += len(good)
    if not rows:
        issues.append('nessuna esecuzione registrata')
    return issues, paired, unpaired, statuses


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'results/results.csv'
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    issues, paired, unpaired, statuses = assess(rows)
    print(f'Coppie NW/Hirschberg confrontate: {paired}')
    print(f'Esecuzioni valide senza confronto score: {unpaired}')
    print(f'Stati: {dict(statuses)}')
    if issues:
        for message in issues[:20]:
            print('ERRORE:', message)
        return 1
    print('Controlli superati; i casi non confrontati non certificano da soli l’ottimalità.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
