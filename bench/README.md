# Benchmark Needleman–Wunsch / Hirschberg

Versione 2026-10-06. Confronto con scoring match +1, mismatch −1, gap −1.
I file `src/` e `include/` sono quelli originali: gli algoritmi, lo scambio interno
in Hirschberg e il parser FASTA della libreria non sono modificati.

## Installazione e avvio

Dipendenze su Ubuntu/WSL: GCC, Make, Python 3 (>=3.9), GNU time. Nessun pacchetto pip.
Se mancanti:

```bash
sudo apt install build-essential python3 time
```

Copiare questa cartella `bench/` nella radice del progetto, accanto a `src/` e
`include/`. Non occorre cambiare il Makefile principale. Usare il Makefile locale:

```bash
make -C bench
make -C bench test
make -C bench quick
```

Il pacchetto completo contiene anche copie inalterate dei sorgenti originali,
quindi si può estrarre in una nuova directory e usare gli stessi comandi.
Non sostituire le proprie copie di `src/` e `include/` se sono state aggiornate
rispetto all'archivio fornito il 6 ottobre.

`make -C bench` compila `bench/align_bench` con
`-O2 -Wall -Wextra -Wpedantic -std=c99` e registra comando, versione compilatore,
hash dei sorgenti e del binario in `build.json`. Ricompila anche quando si cambiano
i flag. Il programma dimostrativo `src/main.c` non viene collegato.

La compilazione va fatta sul proprio WSL ARM64: nessun eseguibile precompilato è
incluso. Se il Makefile principale contiene i vecchi target, possono restare,
ma i comandi documentati qui usano esclusivamente `make -C bench`.

## FASTA originali e preprocessing

Il pacchetto **non include i FASTA biologici scaricati**. Sono predisposte le quattro
coppie concordate in `datasets.csv`:

| dataset_id | file1 | file2 |
|---|---|---|
| rrna16s | NR_112059.1.fasta | NR_116005.1.fasta |
| hiv1 | NC_001802.1.fasta | U46016.1.fasta |
| sars_cov2 | NC_045512.2.fasta | OK091006.1.fasta |
| plastidi | NC_000932.1.fasta | NC_001320.1.fasta |

Scaricare da NCBI Nucleotide gli accession con versione indicati, in formato FASTA,
un record per file; controllarne l'intestazione. Salvare i file con questi nomi
in `bench/data/raw/`. Non sono inferite lunghezze a partire dagli accession: saranno
lette dai file effettivamente forniti. Per file con nomi diversi, modificare il CSV
oppure rinominarli.

Poi, una sola volta:

```bash
make -C bench prepare
```

Il preprocessing:

- accetta A, C, G, T e N, anche minuscoli, sequenze su più righe e spazi bianchi;
- rifiuta record multipli, sequenze vuote e altre ambiguità IUPAC (R, Y, ecc.);
- lascia invariati i file originali;
- sostituisce esclusivamente le N, scegliendo una base tra A, C, G, T;
- usa xorshift64* con seed principale predefinito `20261006`;
- deriva il seed di ogni file dai primi 8 byte big-endian di
  `SHA256(str(seed_principale) + ':' + sha256_file_originale)`, così l'ordine
  dei dataset non cambia il preprocessing; seed zero usa lo stato non nullo
  `0x9e3779b97f4a7c15`;
- scrive FASTA ACGT in `bench/data/processed/` e `manifest.json` con intestazioni,
  lunghezze, numero di N sostituite, seed e hash SHA-256 originali/preprocessati.

Il runner verifica hash, alfabeto e lunghezze prima delle misure e riusa gli stessi
file in tutte le ripetizioni. I seed di preprocessing sono nel manifest; nella
riga CSV FASTA `seed=0` significa che il processo C non genera input casuali.

Per evitare modifiche accidentali ai dati di una campagna, una directory di output
esistente non viene sovrascritta. Per una seconda preparazione scegliere un'altra:

```bash
python3 bench/prepare_fasta.py --seed 20261006 --out bench/data/processed_v2
python3 bench/run_benchmark.py real --manifest bench/data/processed_v2/manifest.json
```

Non aggiungere BASE_N e non cambiare lo scoring: il preprocessamento è esterno
agli algoritmi. Lo score dei dati reali serve al confronto algoritmico, non a
un'interpretazione biologica.

## Campagne

```bash
make -C bench quick       # verifica funzionale sintetica
make -C bench full        # campagna sintetica completa
make -C bench real        # soltanto i dataset del manifest
make -C bench all-data    # full + real nella stessa campagna
```

`all-data` e `real` richiedono che il preprocessing sia già completato. Tutti gli
input reali vengono controllati prima di iniziare, anche con `all-data`.
Ogni campagna crea una **nuova** sottodirectory datata in `bench/results/`.
Non vengono mai sovrascritti risultati precedenti. Check degli score e riepilogo
statistico sono eseguiti automaticamente alla fine.

| Profilo | Lunghezze quadrate | Famiglie | Istanze/ripetizioni |
|---|---|---|---|
| quick | 250, 500, 1000, 2000 | random, mutated | 3 |
| full | 250, 500, 1000, 2000, 4000, 8000, 12000 | random, mutated, identical | 7 |
| real | lunghezze effettive FASTA | fasta | 7 sulla stessa coppia |
| all | full + real | tutte le precedenti | 7 |

- `random`: due sequenze generate pseudocasualmente.
- `mutated`: seconda sequenza derivata dalla prima con sole sostituzioni,
  probabilità predefinita 0.05; la nuova base differisce dalla vecchia.
- `identical`: copia esatta; richiede n=m, come `mutated`.
- I sintetici usano seed `1000 + rep`, quindi **istanze diverse**, non ripetizioni
  della medesima coppia. NW e Hirschberg condividono sempre l'istanza.
- L'ordine è NW/Hirschberg nelle ripetizioni dispari e Hirschberg/NW nelle pari.
- Warm-up: una chiamata 128×128 per algoritmo in processi separati prima della
  campagna, esclusa dai risultati. Non mantiene calde cache private tra processi.
- Tutte le misure sono sequenziali: nessuna parallelizzazione dei benchmark.

### Input rettangolari

`quick`: 500×2000 e 2000×500.
`full`: 1000×10000, 10000×1000, 2000×20000, 20000×2000.

Si genera la coppia con le dimensioni canoniche (corta, lunga), seed
`SEED_BASE + 100 + rep`; la prova inversa usa **esattamente gli stessi input
scambiati** tramite `--swap`. Non si rigenerano semplicemente input con dimensioni
invertite. I risultati dei due orientamenti rimangono separati.

Il wrapper originale di Hirschberg dispone già lo strand più corto sulle righe:
i vettori DP dipendono da min(n,m). Il precedente README che affermava il contrario
è stato corretto. RSS include anche input/output e overhead, quindi non è una misura
diretta della sola memoria dei vettori DP.

## Misure e memoria

- Wall time: CLOCK_MONOTONIC; CPU time: CLOCK_PROCESS_CPUTIME_ID.
- Il timer circonda la chiamata all'algoritmo: include le sue allocazioni e la
  ricostruzione dell'allineamento. Generazione, lettura FASTA e verifica sono fuori.
- `max_rss_kb`: GNU `/usr/bin/time -q -f %M`, unità **KiB** su Linux/WSL.
  Riguarda l'intero processo C: anche parsing, input, output e validazione.
  Non include il processo Python e non rappresenta RAM + swap.
- I file temporanei di GNU time e stdout sono separati. Gli errori non vengono
  concatenati dentro le righe CSV.
- Piccole esecuzioni possono mostrare misure RSS poco informative, anche zero
  secondo il sistema e la granularità del rilevamento. Non vengono inventati
  valori sostitutivi. I test qui forniti non sono misure prestazionali di tesi.

### Budget Needleman–Wunsch

La soglia predefinita è **7 GiB**, scelta come margine rispetto ai 9.7 GiB visibili
nel WSL indicato dall'utente. Si calcola:

```text
(n+1)*(m+1)*(sizeof(int)+sizeof(Direction)) + 3*(n+m)*sizeof(Base)
```

Le prime due componenti sono le matrici score/pred; l'ultimo termine approssima
input e capacità degli output. Le dimensioni dei tipi vengono lette dal binario
con `--abi` (qui normalmente 4, 1, 1 byte).

**È una stima e una soglia di ammissione, non un limite imposto al processo**:
non include ogni overhead e non garantisce l'assenza di OOM. Quando supera la
soglia, NW non parte e si registra `skipped_memory_budget`; Hirschberg viene
comunque eseguito. Per plastidi di circa 150000×150000 questo evita un tentativo
impraticabile. Lo stato non significa che sia stato osservato un OOM.

Personalizzare, mantenendo traccia del valore in tesi:

```bash
python3 bench/run_benchmark.py all --nw-max-gib 7
```

`--nw-max-gib 0` disattiva il filtro. La soglia non limita Hirschberg.
Il timeout predefinito è disattivato; Hirschberg resta quadratico nel tempo anche
quando la memoria è contenuta. Un timeout opzionale:

```bash
python3 bench/run_benchmark.py real --timeout 3600
```

Si applica al processo completo, incluso parsing/verifica, non alla sola regione
temporizzata. Il timeout termina l'intero gruppo di processi e viene registrato
come `timeout`; non si dichiara un OOM a partire da un generico codice di uscita.

## Risultati e controlli

Ogni directory contiene:

- `results.csv`: misure individuali, inclusi fallimenti e casi saltati;
- `summary.csv`: gruppi per dataset, algoritmo, dimensioni, famiglia, tasso e
  orientamento, con mediana, media e deviazione standard **campionaria**;
- `errors.log`: comandi eseguiti e stderr;
- `metadata.json`: parametri, scoring, piattaforma, memoria/swap osservata,
  compilatore/flag, hash sorgenti/binario/script, manifest FASTA e stato campagna.

`rep` distingue ogni ripetizione. `pair_id` collega i due algoritmi nella stessa
campagna: non concatenare CSV di campagne diverse senza rendere univoci gli ID.
`dataset_id` impedisce di fondere due coppie FASTA diverse ma di uguale lunghezza.
`orientation` distingue forward/reverse. Gli hash FASTA identificano gli input;
per sintetici valgono seed, modalità, dimensioni, orientamento e versione sorgenti.

Il processo C verifica la ricostruzione degli input rimuovendo i gap, le lunghezze
degli allineamenti, l'assenza di colonne gap-gap e il ricalcolo dello score.
Il runner confronta poi gli score NW/Hirschberg **per ciascuna ripetizione**:
l'allineamento può essere diverso a parità di score.

| Stato | Significato | Entra nelle statistiche? |
|---|---|---|
| ok | esecuzione valida | sì |
| skipped_memory_budget | NW non avviato per stima sopra soglia | no |
| invalid_alignment | controllo strutturale/score interno fallito | no |
| score_mismatch | i due algoritmi danno score diversi | no, entrambi esclusi |
| timeout | superata durata massima | no |
| failed_N | codice N dal comando, causa in errors.log | no |
| invalid_output / invalid_rss | protocollo di misura non valido | no |

`score_check=matched` indica un confronto riuscito;
`not_compared` segnala che l'altro algoritmo non ha fornito un risultato valido.
Un allineamento valido di Hirschberg può essere incluso anche quando NW è stato
saltato, ma in quel caso **l'ottimalità non è verificata dal confronto tra algoritmi**.

Il riepilogo conserva gruppi con zero misure valide, senza inventare medie.
Mostra tentativi, misure valide, esclusi, saltati, confronti e conteggi per stato.
La deviazione standard è lasciata vuota quando ci sono meno di due misure valide.
Le statistiche sintetiche includono variabilità degli input e del sistema; quelle
FASTA riguardano ripetizioni del medesimo input.

Il verificatore restituisce errore anche per allineamenti invalidi, fallimenti,
timeout, duplicati o righe mancanti; gli skip deliberati sono ammessi e conteggiati.
In caso di Ctrl+C il runner termina il figlio e conserva le righe già completate;
`metadata.json` segnala l'interruzione. Non c'è ripresa automatica di campagne.

Per rieseguire i controlli su una campagna (sostituire il percorso):

```bash
python3 bench/check_scores.py bench/results/NOME_CAMPAGNA/results.csv
python3 bench/summarize.py bench/results/NOME_CAMPAGNA/results.csv bench/results/NOME_CAMPAGNA/summary.csv
```

## Personalizzazione

```bash
python3 bench/run_benchmark.py full --reps 3 --sizes 500 1000 2000 --modes random mutated
python3 bench/run_benchmark.py quick --no-rectangles --out-dir bench/results/prova_1
```

`--out-dir` deve indicare una directory nuova. Rimangono supportate le variabili
REPS, SIZES, MODES, SEED_BASE, MUTATION_RATE, BENCH, NW_MAX_GIB, TIMEOUT_S.
OUT e ERROR_LOG del vecchio runner sono rifiutate: usare `--out-dir`.
La probabilità di mutazione nel generatore conserva la risoluzione di 10^-6
dell'implementazione iniziale.

Lo script shell rimane disponibile:

```bash
bash bench/run_benchmark.sh quick
```

Per chiamare soltanto il binario C (nessuna raccolta RSS o confronto automatico):

```bash
bench/align_bench nw 1000 1000 42 mutated 0.05
bench/align_bench hirschberg 1000 10000 42 random --swap
bench/align_bench nw --fasta bench/data/processed/seq_01.fasta bench/data/processed/seq_02.fasta
```

I FASTA diretti devono essere già ACGT; una N viene rifiutata dal parser originale.
Il comando C produce dieci colonne senza header; è il runner ad aggiungere RSS,
stato e metadati. Eseguire il binario da solo non sostituisce la campagna.

## Collegamento con il capitolo 4

Il codice realizza il protocollo discusso: 7 istanze sintetiche per configurazione,
7 ripetizioni FASTA, ordine alternato, preprocessing una volta sola, tempo della
chiamata e RSS del processo intero. Nel testo finale riportare anche l'effettiva
soglia NW (default 7 GiB) e l'eventuale timeout. Distinguere i casi non avviati per
stima di memoria dai fallimenti osservati. Nei risultati usare lunghezze effettive
nel manifest, non le approssimazioni degli accession.

L'ultimo PDF fornito in questa sessione precede le modifiche dichiarate al capitolo
4; il pacchetto implementa il protocollo concordato in conversazione e non modifica
il documento della tesi.
