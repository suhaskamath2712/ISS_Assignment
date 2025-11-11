#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include <string>

// Simple CSR (Compressed Sparse Row) matrix structure
struct CSRMatrix {
    int n;                        // Number of rows/columns
    int nnz;                      // Number of non-zero elements
    std::vector<int> row_ptr;     // Row pointers (size n+1)
    std::vector<int> col_idx;     // Column indices (size nnz)
    std::vector<double> values;   // Matrix values (size nnz)
    
    CSRMatrix() : n(0), nnz(0) {}
};

// Function declarations
bool load_matrix(const std::string& filename, CSRMatrix& matrix);
bool load_vector(const std::string& filename, std::vector<double>& vec);
void print_matrix_info(const CSRMatrix& matrix);

#endif // SPARSE_MATRIX_H