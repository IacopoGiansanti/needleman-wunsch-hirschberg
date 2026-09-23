#!/usr/bin/env bash
set -u
set -o pipefail

export LC_ALL=C

BENCH="${BENCH:-./align_bench}"
PROFILE="${1:-full}"
OUT="${OUT:-results/results.csv}"
ERROR_LOG="${ERROR_LOG:-results/errors.log}"
MUTATION_RATE="${MUTATION_RATE:-0.05}"
SEED_BASE="${SEED_BASE:-1000}"

mkdir -p "$(dirname "$OUT")"
: > "$ERROR_LOG"

case "$PROFILE" in
    quick)
        REPS="${REPS:-3}"
        SIZES=(${SIZES:-250 500 1000 2000})
        MODES=(${MODES:-random mutated})
        RECT_PAIRS=("500 2000" "2000 500")
        ;;
    full)
        REPS="${REPS:-7}"
        SIZES=(${SIZES:-250 500 1000 2000 4000 8000 12000})
        MODES=(${MODES:-random mutated identical})
        RECT_PAIRS=(
            "1000 10000"
            "10000 1000"
            "2000 20000"
            "20000 2000"
        )
        ;;
    *)
        echo "Uso: $0 [quick|full]" >&2
        exit 1
        ;;
esac

if [[ ! -x "$BENCH" ]]; then
    echo "Errore: $BENCH non esiste o non e' eseguibile." >&2
    echo "Dal root del progetto esegui prima: make bench" >&2
    exit 1
fi

if [[ ! -x /usr/bin/time ]]; then
    echo "Errore: /usr/bin/time non disponibile." >&2
    echo "Su Ubuntu/WSL: sudo apt install time" >&2
    exit 1
fi

echo "algorithm,n,m,seed,input_mode,mutation_rate,wall_time_s,cpu_time_s,score,valid,max_rss_kb,status" > "$OUT"

run_one() {
    local alg="$1"
    local n="$2"
    local m="$3"
    local seed="$4"
    local mode="$5"

    local tmp_mem tmp_out
    tmp_mem="$(mktemp)"
    tmp_out="$(mktemp)"

    /usr/bin/time -f "%M" -o "$tmp_mem" \
        "$BENCH" "$alg" "$n" "$m" "$seed" "$mode" "$MUTATION_RATE" \
        > "$tmp_out" 2>> "$ERROR_LOG"
    local rc=$?

    if [[ $rc -eq 0 || $rc -eq 2 ]]; then
        local row rss status
        row="$(cat "$tmp_out")"
        rss="$(cat "$tmp_mem")"
        if [[ $rc -eq 0 ]]; then
            status="ok"
        else
            status="invalid_alignment"
        fi
        echo "${row},${rss},${status}" >> "$OUT"
    else
        # In caso di malloc failure / OOM / kill, manteniamo l'istanza nel CSV.
        echo "${alg},${n},${m},${seed},${mode},${MUTATION_RATE},NA,NA,NA,NA,NA,failed_${rc}" >> "$OUT"
    fi

    rm -f "$tmp_mem" "$tmp_out"
}

echo "Profilo: $PROFILE"
echo "Output : $OUT"
echo "Reps   : $REPS"
echo

echo "== Warm-up (non registrato) =="
for alg in nw hirschberg; do
    "$BENCH" "$alg" 128 128 42 random "$MUTATION_RATE" >/dev/null 2>> "$ERROR_LOG" || true
done

echo "== Input quadrati =="
for mode in "${MODES[@]}"; do
    for n in "${SIZES[@]}"; do
        for rep in $(seq 1 "$REPS"); do
            seed=$((SEED_BASE + rep))

            # Alterniamo l'ordine per non favorire sempre lo stesso algoritmo.
            if (( rep % 2 == 1 )); then
                order=(nw hirschberg)
            else
                order=(hirschberg nw)
            fi

            for alg in "${order[@]}"; do
                printf 'mode=%-9s n=%-6s rep=%-2s alg=%s\n' "$mode" "$n" "$rep" "$alg"
                run_one "$alg" "$n" "$n" "$seed" "$mode"
            done
        done
    done
done

echo "== Input rettangolari (random) =="
for pair in "${RECT_PAIRS[@]}"; do
    read -r n m <<< "$pair"

    for rep in $(seq 1 "$REPS"); do
        seed=$((SEED_BASE + 100 + rep))

        if (( rep % 2 == 1 )); then
            order=(nw hirschberg)
        else
            order=(hirschberg nw)
        fi

        for alg in "${order[@]}"; do
            printf 'rect n=%-6s m=%-6s rep=%-2s alg=%s\n' "$n" "$m" "$rep" "$alg"
            run_one "$alg" "$n" "$m" "$seed" random
        done
    done
done

echo
echo "Benchmark completato: $OUT"
echo "Log errori          : $ERROR_LOG"

