#!/usr/bin/env python3
"""Campagna sequenziale: ogni misura proviene da un processo C separato."""
import argparse
import csv
import io
import json
import math
import os
import platform
import signal
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from prepare_fasta import sha256, read_fasta

HERE = Path(__file__).resolve().parent
TIME_TOOL = os.environ.get('GNU_TIME', '/usr/bin/time')
C_FIELDS = ['algorithm', 'n', 'm', 'seed', 'input_mode', 'mutation_rate',
            'wall_time_s', 'cpu_time_s', 'score', 'valid']
FIELDS = C_FIELDS + ['max_rss_kb', 'status', 'dataset_id', 'rep', 'orientation',
                    'pair_id', 'input1_sha256', 'input2_sha256',
                    'nw_estimated_bytes', 'exit_code', 'score_check']


def env(name, default):
    return os.environ.get(name, default)


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('profile', nargs='?', choices=['quick', 'full', 'real', 'all'], default='full')
    p.add_argument('--bench', type=Path, default=env('BENCH', HERE / 'align_bench'))
    p.add_argument('--manifest', type=Path, default=HERE / 'data/processed/manifest.json')
    p.add_argument('--out-dir', type=Path, help='Directory nuova, mai sovrascritta')
    p.add_argument('--reps', type=int, default=env('REPS', None))
    p.add_argument('--sizes', nargs='+', type=int, default=env('SIZES', '').split() or None)
    p.add_argument('--modes', nargs='+', choices=['random', 'mutated', 'identical'],
                   default=env('MODES', '').split() or None)
    p.add_argument('--seed-base', type=int, default=env('SEED_BASE', 1000))
    p.add_argument('--mutation-rate', type=float, default=env('MUTATION_RATE', 0.05))
    p.add_argument('--nw-max-gib', type=float, default=env('NW_MAX_GIB', 7),
                   help='Budget stimato NW (GiB), default 7; 0 disabilita il filtro')
    p.add_argument('--timeout', type=float, default=env('TIMEOUT_S', 0),
                   help='Secondi per processo; 0 = nessun timeout')
    p.add_argument('--no-rectangles', action='store_true')
    a = p.parse_args()
    quick = a.profile == 'quick'
    a.reps = a.reps if a.reps is not None else (3 if quick else 7)
    a.sizes = list(map(int, a.sizes or ([250, 500, 1000, 2000] if quick else
                                      [250, 500, 1000, 2000, 4000, 8000, 12000])))
    a.modes = a.modes or (['random', 'mutated'] if quick else ['random', 'mutated', 'identical'])
    if a.reps < 1 or any(n < 1 or n > 2**30 - 1 for n in a.sizes):
        p.error('reps e dimensioni devono essere positive; dimensioni <= 2^30-1')
    if len(set(a.sizes)) != len(a.sizes) or len(set(a.modes)) != len(a.modes):
        p.error('dimensioni o modalità duplicate')
    if any(m not in ('random', 'mutated', 'identical') for m in a.modes):
        p.error('modalità non valida')
    if not 0 <= a.seed_base < 2**64 - a.reps - 101:
        p.error('seed-base fuori intervallo')
    if not math.isfinite(a.mutation_rate) or not 0 <= a.mutation_rate <= 1:
        p.error('mutation-rate deve appartenere a [0,1]')
    if any(not math.isfinite(x) or x < 0 for x in (a.timeout, a.nw_max_gib)):
        p.error('timeout e budget devono essere finiti e non negativi')
    # Gli script vecchi scrivevano sempre lo stesso CSV: ora nessuna sovrascrittura.
    if 'OUT' in os.environ or 'ERROR_LOG' in os.environ:
        p.error('OUT/ERROR_LOG non più supportati: usare --out-dir NUOVA_DIRECTORY')
    return a


def load_real(path):
    data = json.loads(path.read_text())
    if data.get('schema_version') != 1 or not data.get('datasets'):
        raise ValueError('manifest non supportato o vuoto')
    ids = set()
    for d in data['datasets']:
        if not d['dataset_id'] or d['dataset_id'] in ids:
            raise ValueError('dataset_id mancante o duplicato nel manifest')
        ids.add(d['dataset_id'])
        for side in ('first', 'second'):
            record = d[side]
            source = (path.parent / record['path']).resolve()
            if sha256(source) != record['processed_sha256']:
                raise ValueError(f'Hash FASTA non corrispondente: {source}')
            _, seq = read_fasta(source)
            if 'N' in seq or len(seq) != record['length']:
                raise ValueError(f'FASTA non preprocessato o lunghezza errata: {source}')
            record['resolved_path'] = str(source)
    return data


def cases(a, real):
    if a.profile in ('quick', 'full', 'all'):
        for mode in a.modes:
            for n in a.sizes:
                for rep in range(1, a.reps + 1):
                    yield dict(dataset_id=f'synthetic_{mode}', rep=rep, orientation='forward',
                        n=n, m=n, seed=a.seed_base + rep, input_mode=mode,
                        mutation_rate=a.mutation_rate if mode == 'mutated' else 0,
                        source_n=n, source_m=n)
        rectangles = [(500, 2000)] if a.profile == 'quick' else [(1000, 10000), (2000, 20000)]
        if not a.no_rectangles:
            for n, m in rectangles:
                for rep in range(1, a.reps + 1):
                    for reverse in (False, True):
                        yield dict(dataset_id='synthetic_rectangle', rep=rep,
                            orientation='reverse' if reverse else 'forward',
                            n=m if reverse else n, m=n if reverse else m,
                            source_n=n, source_m=m, seed=a.seed_base + 100 + rep,
                            input_mode='random', mutation_rate=0)
    if real:
        for d in real['datasets']:
            for rep in range(1, a.reps + 1):
                yield dict(dataset_id=d['dataset_id'], rep=rep, orientation='forward',
                    n=d['first']['length'], m=d['second']['length'], seed=0,
                    input_mode='fasta', mutation_rate=0,
                    file1=d['first']['resolved_path'], file2=d['second']['resolved_path'],
                    input1_sha256=d['first']['processed_sha256'],
                    input2_sha256=d['second']['processed_sha256'])


def command(binary, algorithm, case):
    if case['input_mode'] == 'fasta':
        return [str(binary), algorithm, '--fasta', case['file1'], case['file2']]
    cmd = [str(binary), algorithm, str(case['source_n']), str(case['source_m']),
           str(case['seed']), case['input_mode'], str(case['mutation_rate'])]
    if case['orientation'] == 'reverse':
        cmd.append('--swap')
    return cmd


def measure(binary, algorithm, case, abi, budget, timeout, log):
    estimate = ((case['n'] + 1) * (case['m'] + 1) *
                (abi['int_bytes'] + abi['direction_bytes']) +
                3 * (case['n'] + case['m']) * abi['base_bytes'])
    row = {key: case.get(key, '') for key in FIELDS}
    row['mutation_rate'] = f'{float(case["mutation_rate"]):.6f}'
    row.update(algorithm=algorithm, nw_estimated_bytes=estimate,
               status='pending', score_check='not_compared')
    if algorithm == 'nw' and budget and estimate > budget:
        row['status'] = 'skipped_memory_budget'
        return row
    cmd = command(binary, algorithm, case)
    log.write(f'\n{case["pair_id"]} {algorithm}: {json.dumps(cmd)}\n')
    log.flush()
    with tempfile.TemporaryDirectory(prefix='align-rss-') as tmp:
        mem = Path(tmp) / 'rss.txt'
        proc = subprocess.Popen([TIME_TOOL, '-q', '-f', '%M', '-o', str(mem)] + cmd,
                                stdout=subprocess.PIPE, stderr=log, text=True,
                                start_new_session=True, env={**os.environ, 'LC_ALL': 'C'})
        timed_out = False
        try:
            output, _ = proc.communicate(timeout=timeout or None)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            output, _ = proc.communicate()
        except BaseException:
            # Non lasciare un allineamento attivo dopo Ctrl+C.
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.communicate()
            raise
        row['exit_code'] = proc.returncode
        rss = mem.read_text().strip() if mem.exists() else ''
        row['max_rss_kb'] = rss if rss.isdigit() else ''
    if timed_out:
        row['status'] = 'timeout'
        return row
    if proc.returncode not in (0, 2):
        row['status'] = f'failed_{proc.returncode}'
        return row
    try:
        lines = list(csv.reader(io.StringIO(output)))
        if len(lines) != 1 or len(lines[0]) != len(C_FIELDS):
            raise ValueError('CSV malformato')
        result = dict(zip(C_FIELDS, lines[0]))
        for key in ('algorithm', 'input_mode'):
            if result[key] != str(row[key]):
                raise ValueError(f'{key} inatteso')
        for key in ('n', 'm', 'seed'):
            if int(result[key]) != int(row[key]):
                raise ValueError(f'{key} inatteso')
        if abs(float(result['mutation_rate']) - float(row['mutation_rate'])) > 0.00000051:
            raise ValueError('mutation_rate inatteso')
        for key in ('wall_time_s', 'cpu_time_s'):
            if not math.isfinite(float(result[key])) or float(result[key]) < 0:
                raise ValueError('tempo non valido')
        int(result['score'])
        if result['valid'] not in ('0', '1'):
            raise ValueError('campo valid non valido')
        row.update(result)
        row['status'] = 'ok' if result['valid'] == '1' and proc.returncode == 0 else 'invalid_alignment'
        if row['status'] == 'ok' and not row['max_rss_kb']:
            row['status'] = 'invalid_rss'
    except (ValueError, KeyError) as e:
        log.write(f'Output non valido: {e}; {output!r}\n')
        row['status'] = 'invalid_output'
    return row


def file_text(path):
    try:
        return Path(path).read_text()
    except OSError:
        return None


def main():
    a = arguments()
    a.bench = a.bench.resolve()
    if not a.bench.is_file() or not os.access(a.bench, os.X_OK):
        raise ValueError('Eseguibile assente: eseguire make -C bench')
    if not Path(TIME_TOOL).is_file():
        raise ValueError('/usr/bin/time assente (su Ubuntu: sudo apt install time)')
    abi = json.loads(subprocess.check_output([str(a.bench), '--abi'], text=True))
    real = load_real(a.manifest.resolve()) if a.profile in ('real', 'all') else None
    build_path = a.bench.parent / 'build.json'
    build = json.loads(build_path.read_text()) if build_path.exists() else None
    binary_hash = sha256(a.bench)
    if build is None or build['binary_sha256'] != binary_hash:
        raise ValueError('build.json assente o incoerente: ricompilare con make -C bench')
    for name, digest in build['sources'].items():
        if sha256(a.bench.parent / name) != digest:
            raise ValueError(f'Sorgente modificato dopo la compilazione: {name}; ricompilare')
    planned = list(cases(a, real))
    for index, case in enumerate(planned, 1):
        case['pair_id'] = f'pair_{index:05d}'
        if case['n'] + case['m'] > 2**31 - 1:
            raise ValueError('Input troppo grande per gli score int')
    out = a.out_dir or HERE / 'results' / (datetime.now().strftime('%Y%m%d_%H%M%S_%f_') + a.profile)
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    meta = dict(started_at=datetime.now(timezone.utc).isoformat(), status='running',
        parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
        platform=platform.platform(), machine=platform.machine(), python=sys.version,
        cpuinfo=file_text('/proc/cpuinfo'), meminfo=file_text('/proc/meminfo'),
        swap=file_text('/proc/swaps'), build=build, abi=abi, scoring=dict(match=1,mismatch=-1,gap=-1),
        gnu_time=subprocess.check_output([TIME_TOOL, '--version'], text=True),
        scripts_sha256={p.name: sha256(p) for p in sorted(HERE.glob('*.py'))},
        planned_pairs=len(planned), real_manifest=real)
    metadata = out / 'metadata.json'
    metadata.write_text(json.dumps(meta, indent=2) + '\n')
    print(f'Profilo: {a.profile}; istanze/ripetizioni: {a.reps}; output: {out}', flush=True)
    try:
        with (out / 'errors.log').open('w') as log, (out / 'results.csv').open('w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=FIELDS)
            writer.writeheader()
            csvfile.flush()
            for alg in ('nw', 'hirschberg'):
                subprocess.run([str(a.bench), alg, '128', '128', '42', 'random'],
                               stdout=subprocess.DEVNULL, stderr=log, check=True, timeout=30)
            for index, case in enumerate(planned, 1):
                pair = []
                try:
                    order = ('nw', 'hirschberg') if case['rep'] % 2 else ('hirschberg', 'nw')
                    for alg in order:
                        print(f'[{index}/{len(planned)}] {case["dataset_id"]} '
                              f'{case["n"]}x{case["m"]} rep={case["rep"]} {alg}', flush=True)
                        pair.append(measure(a.bench, alg, case, abi,
                            a.nw_max_gib * 1024**3, a.timeout, log))
                        print(f'  {pair[-1]["status"]}', flush=True)
                    if all(r['status'] == 'ok' for r in pair):
                        equal = pair[0]['score'] == pair[1]['score']
                        for r in pair:
                            r['score_check'] = 'matched' if equal else 'mismatch'
                            if not equal:
                                r['status'] = 'score_mismatch'
                finally:
                    # Anche se interrotti, conservare il primo risultato completato.
                    writer.writerows(pair)
                    csvfile.flush()
        check = subprocess.run([sys.executable, str(HERE / 'check_scores.py'), str(out / 'results.csv')])
        subprocess.run([sys.executable, str(HERE / 'summarize.py'), str(out / 'results.csv'),
                        str(out / 'summary.csv')], check=True)
        meta['status'] = 'completed' if check.returncode == 0 else 'completed_with_errors'
        return check.returncode
    except BaseException:
        meta['status'] = 'interrupted_or_failed'
        raise
    finally:
        meta['finished_at'] = datetime.now(timezone.utc).isoformat()
        metadata.write_text(json.dumps(meta, indent=2) + '\n')
        print(f'Risultati salvati in {out}', flush=True)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(130)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as e:
        sys.exit(f'Errore: {e}')
