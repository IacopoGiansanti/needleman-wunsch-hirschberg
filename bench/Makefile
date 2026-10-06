CC = gcc
CFLAGS = -O2 -Wall -Wextra -Wpedantic -std=c99
CPPFLAGS = -I../include
PYTHON = python3
SOURCES = benchmark.c bench_adapter.c ../src/alignment.c ../src/base.c ../src/fasta.c ../src/strand.c ../src/scoring.c ../src/middle_edge.c ../src/Needleman-Wunsch.c ../src/Hirschberg.c
HEADERS = $(wildcard ../include/*.h) bench_adapter.h

.PHONY: all bench quick full real all-data prepare test clean
all: bench
# Ricompila anche quando cambiano CC/CFLAGS, salvando i flag effettivi.
bench:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) $(LDFLAGS) $(LDLIBS) -o align_bench
	$(PYTHON) record_build.py --compiler "$(CC)" --command '$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) $(LDFLAGS) $(LDLIBS) -o align_bench' $(SOURCES) $(HEADERS)
quick: bench
	$(PYTHON) run_benchmark.py quick
full: bench
	$(PYTHON) run_benchmark.py full
real: bench
	$(PYTHON) run_benchmark.py real
all-data: bench
	$(PYTHON) run_benchmark.py all
prepare:
	$(PYTHON) prepare_fasta.py
test: bench
	$(PYTHON) -m unittest discover -s tests -v
clean:
	rm -f align_bench build.json
