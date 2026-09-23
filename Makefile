CC = gcc

CFLAGS = -O2 -Wall -Wextra -Wpedantic -std=c99
CPPFLAGS = -Iinclude

SRC = src/base.c \
	  src/strand.c \
	  src/alignment.c \
	  src/scoring.c \
	  src/fasta.c \
	  src/middle_edge.c \
	  src/Needleman-Wunsch.c \
	  src/Hirschberg.c \
	  src/main.c

OBJ = $(SRC:.c=.o)

TARGET = align


# =========================
# Main program
# =========================

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


# =========================
# Benchmark
# =========================

BENCH_TARGET = bench/align_bench

BENCH_CFLAGS = -O3 -DNDEBUG -Wall -Wextra -Wpedantic -std=c99

CORE_SRC = src/base.c \
	   src/strand.c \
	   src/alignment.c \
	   src/scoring.c \
	   src/fasta.c \
	   src/middle_edge.c \
	   src/Needleman-Wunsch.c \
	   src/Hirschberg.c

BENCH_SRC = bench/benchmark.c \
	    bench/bench_adapter.c \
	    $(CORE_SRC)

bench: $(BENCH_TARGET)

$(BENCH_TARGET): $(BENCH_SRC)
	$(CC) $(BENCH_CFLAGS) $(CPPFLAGS) -Ibench $(BENCH_SRC) -o $(BENCH_TARGET)

benchmark-quick: bench
	cd bench && bash run_benchmark.sh quick

benchmark-full: bench
	cd bench && bash run_benchmark.sh full

benchmark-check:
	cd bench && python3 check_scores.py results/results.csv

benchmark-summary:
	cd bench && python3 summarize.py results/results.csv results/summary.csv


# =========================
# Cleaning
# =========================

clean:
	rm -f $(OBJ) $(TARGET) $(BENCH_TARGET)


.PHONY: bench benchmark-quick benchmark-full benchmark-check benchmark-summary clean
