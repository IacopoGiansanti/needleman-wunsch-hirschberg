CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -Wpedantic -std=c99
CPPFLAGS ?= -I../include

CORE_SRC = ../src/base.c \
           ../src/strand.c \
           ../src/alignment.c \
           ../src/scoring.c \
           ../src/fasta.c \
           ../src/middle_edge.c \
           ../src/Needleman-Wunsch.c \
           ../src/Hirschberg.c

TARGET = alignment_dump

all: $(TARGET)

$(TARGET): alignment_dump.c $(CORE_SRC)
	$(CC) $(CFLAGS) $(CPPFLAGS) alignment_dump.c $(CORE_SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all clean
