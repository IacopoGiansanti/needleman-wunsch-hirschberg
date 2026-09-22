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

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: clean
