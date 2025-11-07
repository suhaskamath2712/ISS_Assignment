	## GNU Makefile to build sequential and OpenMP variants of k-means (Linux)

# Compiler and common flags
CC           ?= gcc
CSTD         ?= c11
BASE_CFLAGS  ?= -O2 -Wall -Wextra -pedantic -std=$(CSTD)
LDLIBS       ?= -lm
OMP_FLAGS    ?= -fopenmp

# Programs and sources
SEQ_BIN := kmeans_sequential
PAR_BIN := kmeans_parallel

SEQ_SRC := k_means_sequential.c
PAR_SRC := k_means_parallel.c

SEQ_OBJ := $(SEQ_SRC:.c=.o)
PAR_OBJ := $(PAR_SRC:.c=.o)

# Build both by default
all: $(SEQ_BIN) $(PAR_BIN)

# Debug build: rebuild with -g
debug: BASE_CFLAGS := -g -Wall -Wextra -pedantic -std=$(CSTD)
debug: clean all

# Sequential compile and link (no OpenMP flags)
$(SEQ_OBJ): $(SEQ_SRC)
	$(CC) $(BASE_CFLAGS) -c $< -o $@

$(SEQ_BIN): $(SEQ_OBJ)
	$(CC) $(BASE_CFLAGS) $^ -o $@ $(LDLIBS)

# Parallel (OpenMP) compile and link
$(PAR_OBJ): $(PAR_SRC)
	$(CC) $(BASE_CFLAGS) $(OMP_FLAGS) -c $< -o $@

$(PAR_BIN): $(PAR_OBJ)
	$(CC) $(BASE_CFLAGS) $(OMP_FLAGS) $^ -o $@ $(LDLIBS) $(OMP_FLAGS)

# Clean artifacts
clean:
	rm -f $(SEQ_OBJ) $(PAR_OBJ) $(SEQ_BIN) $(PAR_BIN)
