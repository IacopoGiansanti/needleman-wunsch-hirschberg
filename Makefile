CC = gcc
CFLAGS = -O2 -Wall -Wextra -Wpedantic -std=c99

SRC = base.c \
	  strand.c \
	  alignment.c \
	  scoring.c \
	  fasta.c \
	  middle_edge.c \
	  Needleman-Wunsch.c \
      Hirschberg.c \
	  main.c

OBJ = $(SRC:.c=.o)
TARGET = align

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

BENCH_TARGET = bench/align_bench
BENCH_CFLAGS = -O3 -DNDEBUG -Wall -Wextra -Wpedantic -std=c99

CORE_SRC = base.c \
	strand.c \
	alignment.c \
	scoring.c \
	fasta.c \
	middle_edge.c \
	Needleman-Wunsch.c \
	Hirschberg.c

BENCH_SRC = bench/benchmark.c \
	bench/bench_adapter.c \
	$(CORE_SRC)

bench: $(BENCH_TARGET)

$(BENCH_TARGET): $(BENCH_SRC)
	$(CC) $(BENCH_CFLAGS) -I. -Ibench $(BENCH_SRC) -o $(BENCH_TARGET)

benchmark-quick: bench
	cd bench && bash run_benchmark.sh quick

benchmark-full: bench
	cd bench && bash run_benchmark.sh full

benchmark-check:
	cd bench && python3 check_scores.py results/results.csv

benchmark-summary:
	cd bench && python3 summarize.py results/results.csv results/summary.csv

.PHONY: bench benchmark-quick benchmark-full benchmark-check benchmark-summary
 
clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: clean
