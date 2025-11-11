#include "sparse_matrix.h"
#include <iostream>
#include <vector>
#include <chrono>

using namespace std;

// Simple serial sparse matrix-vector multiplication
void spmv_serial(const CSRMatrix& A, const vector<double>& x, vector<double>& y)
{
    auto start = chrono::high_resolution_clock::now();
    
    for (int i = 0; i < A.n; i++) {
        y[i] = 0.0;
        for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; j++) {
            y[i] += A.values[j] * x[A.col_idx[j]];
        }
    }
    
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
}

int main()
{
    cout << "=== NLPKKT240 Sparse Matrix-Vector Multiplication Example ===" << endl;
    
    // File paths - read from shared directory (not copied locally)
    const string matrix_file = "/scratch/public/sparse_matrix_data/nlpkkt240_matrix.bin";
    const string vector_file = "/scratch/public/sparse_matrix_data/nlpkkt240_vector.bin";
    
    // Load matrix and vector
    CSRMatrix A;
    vector<double> x, y;
    
    cout << "\nLoading matrix from: " << matrix_file << endl;
    if (!load_matrix(matrix_file, A)) {
        cerr << "Failed to load matrix." << endl;
        cerr << "Make sure you have read access to /scratch/public/sparse_matrix_data/" << endl;
        return -1;
    }
    
    cout << "Loading vector from: " << vector_file << endl;
    if (!load_vector(vector_file, x)) {
        cerr << "Failed to load vector." << endl;
        return -1;
    }
    
    // Print matrix information
    print_matrix_info(A);
    
    // Prepare result vector
    y.resize(A.n, 0.0);
    
    // Perform serial SpMV
    cout << "\nPerforming serial sparse matrix-vector multiplication..." << endl;
    spmv_serial(A, x, y);
    
    // Verify result (check a few values)
    cout << "\nSample results:" << endl;
    cout << "y[0] = " << y[0] << endl;
    cout << "y[" << A.n/2 << "] = " << y[A.n/2] << endl;
    cout << "y[" << A.n-1 << "] = " << y[A.n-1] << endl;
    
    return 0;
}
