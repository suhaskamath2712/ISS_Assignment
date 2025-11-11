CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall

# For OpenMP (students will uncomment when implementing)
# CXXFLAGS += -fopenmp

# For MPI (students will modify when implementing MPI version)
CXX = mpicxx

TARGET = main
SOURCES = main.cpp sparse_matrix.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

clean:
	rm -f $(TARGET)

test: $(TARGET)
	./$(TARGET)

help:
	@echo "Build: make"
	@echo "Run:   make test"
	@echo "Clean: make clean"

.PHONY: all clean test help