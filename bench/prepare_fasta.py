#!/usr/bin/env python3
"""Preprocessing offline, una sola volta; manifest con hash e conteggi delle ambiguita IUPAC."""
import argparse
import csv
from collections import Counter
import hashlib
import json
from pathlib import Path

MASK = (1 << 64) - 1
# Codici IUPAC del DNA, ordinati A,C,G,T per una scelta riproducibile.
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_fasta(path, allow_ambiguity=False):
    header = None
    chunks = []
    with open(path, encoding='ascii') as f:
        for number, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    raise ValueError(f'{path}:{number}: più di un record FASTA')
                header = line[1:].strip()
                if not header:
                    raise ValueError(f'{path}:{number}: intestazione vuota')
            else:
                if header is None:
                    raise ValueError(f'{path}:{number}: manca intestazione FASTA')
                seq = ''.join(line.split()).upper()
                bad = set(seq) - (set(IUPAC) if allow_ambiguity else set('ACGT'))
                if bad:
                    raise ValueError(f'{path}:{number}: simboli non supportati: {sorted(bad)}')
                chunks.append(seq)
    seq = ''.join(chunks)
    if header is None or not seq:
        raise ValueError(f'{path}: sequenza vuota o intestazione assente')
    return header, seq


def resolve_ambiguities(seq, seed):
    state = seed or 0x9e3779b97f4a7c15
    out = []
    for base in seq:
        allowed = IUPAC[base]
        if len(allowed) > 1:
            # Rejection sampling: evita bias modulo per insiemi di tre basi.
            limit = (1 << 64) - ((1 << 64) % len(allowed))
            while True:
                state ^= state >> 12
                state ^= (state << 25) & MASK
                state ^= state >> 27
                value = (state * 2685821657736338717) & MASK
                if value < limit:
                    break
            base = allowed[value % len(allowed)]
        out.append(base)
    return ''.join(out)


def prepare(table, raw_dir, output, seed):
    if not 0 <= seed <= MASK:
        raise ValueError('seed fuori dall’intervallo uint64')
    with open(table, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    if not rows or not {'dataset_id', 'file1', 'file2'} <= rows[0].keys():
        raise ValueError('CSV richiesto: dataset_id,file1,file2')
    ids = [r['dataset_id'] for r in rows]
    if any(not x.strip() for x in ids) or len(set(ids)) != len(ids):
        raise ValueError('dataset_id vuoti o duplicati')
    # Validare TUTTI gli input prima di creare la directory di output.
    loaded = {}
    for row in rows:
        for key in ('file1', 'file2'):
            path = (raw_dir / row[key]).resolve()
            if path not in loaded:
                header, seq = read_fasta(path, allow_ambiguity=True)
                source_hash = sha256(path)
                # Seed per contenuto: indipendente da ordine e nome del file.
                derived = int.from_bytes(hashlib.sha256(
                    f'{seed}:{source_hash}'.encode()).digest()[:8], 'big')
                loaded[path] = (header, seq, source_hash, derived)
    output.mkdir(parents=True, exist_ok=False)
    records = {}
    for i, (path, (header, seq, source_hash, derived)) in enumerate(loaded.items(), 1):
        processed = output / f'seq_{i:02d}.fasta'
        clean = resolve_ambiguities(seq, derived)
        counts = {symbol: count for symbol, count in sorted(Counter(seq).items())
                  if symbol not in 'ACGT'}
        with processed.open('w', encoding='ascii', newline='\n') as f:
            f.write(f'>{header}\n')
            for start in range(0, len(clean), 80):
                f.write(clean[start:start + 80] + '\n')
        records[path] = dict(path=processed.name, source_path=str(path),
            source_sha256=source_hash, processed_sha256=sha256(processed),
            length=len(clean), n_replaced=seq.count('N'),
            ambiguity_counts=counts, ambiguity_replaced=sum(counts.values()),
            seed=derived, header=header)
    datasets = []
    for row in rows:
        datasets.append(dict(dataset_id=row['dataset_id'],
            first=records[(raw_dir / row['file1']).resolve()],
            second=records[(raw_dir / row['file2']).resolve()]))
    manifest = dict(schema_version=1, preprocessing='xorshift64star-IUPAC-v2',
                    master_seed=seed, datasets=datasets)
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Preparati {len(datasets)} dataset: {output / "manifest.json"}')
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve().parent
    p.add_argument('--pairs', type=Path, default=here / 'datasets.csv')
    p.add_argument('--raw-dir', type=Path, default=here / 'data/raw')
    p.add_argument('--out', type=Path, default=here / 'data/processed')
    p.add_argument('--seed', type=int, default=20261006)
    a = p.parse_args()
    try:
        prepare(a.pairs, a.raw_dir, a.out, a.seed)
    except (ValueError, OSError, UnicodeError) as e:
        p.exit(1, f'Errore: {e}\n')


if __name__ == '__main__':
    main()
