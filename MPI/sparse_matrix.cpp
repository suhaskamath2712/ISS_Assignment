#include "sparse_matrix.h"
#include <iostream>
#include <fstream>

// Load CSR matrix from binary file
bool load_matrix(const std::string& filename, CSRMatrix& matrix) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open matrix file " << filename << std::endl;
        return false;
    }
    
    // Read dimensions
    file.read(reinterpret_cast<char*>(&matrix.n), sizeof(int));
    file.read(reinterpret_cast<char*>(&matrix.nnz), sizeof(int));
    
    // Resize vectors
    matrix.row_ptr.resize(matrix.n + 1);
    matrix.col_idx.resize(matrix.nnz);
    matrix.values.resize(matrix.nnz);
    
    // Read data
    file.read(reinterpret_cast<char*>(matrix.row_ptr.data()), (matrix.n + 1) * sizeof(int));
    file.read(reinterpret_cast<char*>(matrix.col_idx.data()), matrix.nnz * sizeof(int));
    file.read(reinterpret_cast<char*>(matrix.values.data()), matrix.nnz * sizeof(double));
    
    file.close();
    
    //std::cout << "Matrix loaded: " << matrix.n << "x" << matrix.n 
    //          << ", nnz = " << matrix.nnz << std::endl;
    return true;
}

// Load vector from binary file
bool load_vector(const std::string& filename, std::vector<double>& vec) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open vector file " << filename << std::endl;
        return false;
    }
    
    // Read vector size
    int n;
    file.read(reinterpret_cast<char*>(&n), sizeof(int));
    
    // Resize and read vector data
    vec.resize(n);
    file.read(reinterpret_cast<char*>(vec.data()), n * sizeof(double));
    
    file.close();
    
    //std::cout << "Vector loaded: size = " << n << std::endl;
    return true;
}

// Print basic matrix information
void print_matrix_info(const CSRMatrix& matrix) {
    std::cout << "\n=== Matrix Information ===" << std::endl;
    std::cout << "Dimensions: " << matrix.n << " x " << matrix.n << std::endl;
    std::cout << "Non-zeros: " << matrix.nnz << std::endl;
    std::cout << "Sparsity: " << (1.0 - (double)matrix.nnz / ((double)matrix.n * matrix.n)) * 100 << "%" << std::endl;
    std::cout << "Memory usage: ~" << 
        (matrix.row_ptr.size() * sizeof(int) + 
         matrix.col_idx.size() * sizeof(int) + 
         matrix.values.size() * sizeof(double)) / (1024.0 * 1024.0) << " MB" << std::endl;
}