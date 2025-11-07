	## GNU Makefile to build sequential and OpenMP variants of k-means (Linux)

# Compiler and common flags
CC           ?= gcc
CSTD         ?= c11
BASE_CFLAGS  ?= -O2 -Wall -Wextra -pedantic -std=$(CSTD)
LDLIBS       ?= -lm
OMP_FLAGS    ?= -fopenmp

# Programs and sources
PAR_BIN := kmeans_parallel
PAR_SRC := k_means_parallel.c
PAR_OBJ := $(PAR_SRC:.c=.o)

# Build both by default
all: $(PAR_BIN)

# Debug build: rebuild with -g
debug: BASE_CFLAGS := -g -Wall -Wextra -pedantic -std=$(CSTD)
debug: clean all

# Parallel (OpenMP) compile and link
$(PAR_OBJ): $(PAR_SRC)
	$(CC) $(BASE_CFLAGS) $(OMP_FLAGS) -c $< -o $@

$(PAR_BIN): $(PAR_OBJ)
	$(CC) $(BASE_CFLAGS) $(OMP_FLAGS) $^ -o $@ $(LDLIBS) $(OMP_FLAGS)

# Clean artifacts
clean:
	rm -f $(PAR_OBJ) $(PAR_BIN)
