# Benchmark Needleman-Wunsch vs Hirschberg

Versione adattata alle API reali del progetto.

Le funzioni del benchmark lavorano direttamente con `Strand` e con `scoringModel`.

# Benchmark Needleman-Wunsch vs Hirschberg

Questa cartella contiene il benchmark separato dal programma dimostrativo `main.c`.
La regione temporizzata contiene **solo** la chiamata all'algoritmo; generazione degli
input e verifica dell'allineamento sono escluse dal tempo misurato.

## Compilazione

Dal root del progetto:

```bash
make bench
```

Viene creato:

```text
bench/align_bench
```

Il benchmark usa lo stesso scoring model del programma principale:

- match = +1
- mismatch = -1
- gap = -1

## Test rapido

Prima di una campagna lunga:

```bash
make benchmark-quick
```

oppure:

```bash
cd bench
./run_benchmark.sh quick
```

## Benchmark completo

```bash
make benchmark-full
```

Il runner esegue piu' ripetizioni per ogni configurazione, alternando l'ordine
Needleman-Wunsch / Hirschberg e usando gli stessi seed per entrambi.

Produce:

```text
bench/results/results.csv
bench/results/errors.log
```

Le colonne principali sono:

- `wall_time_s`: tempo reale della sola chiamata all'algoritmo;
- `cpu_time_s`: CPU time della sola chiamata all'algoritmo;
- `max_rss_kb`: picco RSS dell'intero processo, misurato con `/usr/bin/time`;
- `score`: score ottimo restituito;
- `valid`: validita' strutturale dell'allineamento restituito;
- `status`: `ok`, `invalid_alignment` oppure un codice di fallimento/OOM.

| Colonna         | Significato                                        |
| --------------- | -------------------------------------------------- |
| `algorithm`     | `nw` oppure `hirschberg`                           |
| `n`             | lunghezza della prima sequenza                     |
| `m`             | lunghezza della seconda sequenza                   |
| `seed`          | seed usato per generare le sequenze casuali        |
| `input_mode`    | `random`, `mutated`, ecc.                          |
| `mutation_rate` | probabilità di mutazione usata per input `mutated` |
| `wall_time_s`   | tempo reale trascorso, in secondi                  |
| `cpu_time_s`    | tempo CPU consumato dal processo                   |
| `score`         | score dell'allineamento ottimo                     |
| `valid`         | `1` se l'allineamento prodotto è valido            |
| `max_rss_kb`    | massimo uso di memoria residente, in KiB           |
| `status`        | `ok` se l'esecuzione è terminata correttamente     |


## Controllo di correttezza

```bash
make benchmark-check
```

Per la stessa istanza, NW e Hirschberg devono ottenere lo stesso score. Non e'
richiesto che producano lo stesso identico alignment, perche' possono esistere piu'
cammini ottimi in presenza di pareggi.

## Riassunto statistico

```bash
make benchmark-summary
```

Crea:

```text
bench/results/summary.csv
```

con mediana, media e deviazione standard di tempo e memoria.

## Famiglie di input

- `random`: sequenze DNA indipendenti;
- `mutated`: la seconda sequenza e' una copia della prima con il 5% di sostituzioni;
- `identical`: le due sequenze coincidono.

## Input rettangolari

Il benchmark completo include anche coppie `n x m` e `m x n`.
Questo e' particolarmente utile per **questa implementazione** di Hirschberg: le
funzioni `nw_score` e `nw_score_reverse` allocano vettori di lunghezza `n + 1`, dove
`n` e' la lunghezza del primo strand. Quindi la memoria lineare e' rispetto al primo
argomento e non viene automaticamente minimizzata scambiando le due sequenze.
Il benchmark rettangolare rende questa caratteristica direttamente osservabile.

## Personalizzazione

Esempi:

```bash
cd bench
REPS=3 SIZES="500 1000 2000 4000" ./run_benchmark.sh full
```

```bash
MODES="random mutated" ./run_benchmark.sh full
```

```bash
OUT=results/esperimento_1.csv ./run_benchmark.sh full
```

## Nota sulla memoria

`max_rss_kb` e' una misura reale del picco di memoria residente del processo, quindi
include anche runtime, input e output. E' la misura principale da usare per il confronto
sperimentale. Per input piccoli il costo fisso del processo puo' mascherare parte della
crescita lineare di Hirschberg; il fenomeno diventa molto piu' chiaro aumentando `n`.

