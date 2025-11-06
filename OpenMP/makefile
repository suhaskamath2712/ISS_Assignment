# GNU Makefile for k_means_clustering.c (Linux)

# Compiler and flags
CC      ?= gcc
CSTD    ?= c11
CFLAGS  ?= -O2 -Wall -Wextra -pedantic -std=$(CSTD)
LDFLAGS ?=
LDLIBS  ?= -lm

# Sources and targets
SRC   := k_means_clustering.c
OBJ   := $(SRC:.c=.o)
BIN   := kmeans


.PHONY: all debug clean run

# Default build: optimized
all: $(BIN)

# Debug build: add symbols and disable optimizations
debug: CFLAGS := -g -O0 -Wall -Wextra -pedantic -std=$(CSTD)
debug: clean $(BIN)

# Link
$(BIN): $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# Compile
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
