/*
================================================================================
Sparse Matrix-Vector Multiplication (SpMV) with MPI Parallelization
================================================================================
This version outputs the results in CSV format:
number_of_processes,computation_time_ms
================================================================================
*/

#include "sparse_matrix.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <fstream> // Needed for file I/O in parallel
#include <mpi.h>

using namespace std;

// File paths for matrix and vector (must be accessible to all ranks)
const string matrix_file = "/scratch/public/sparse_matrix_data/nlpkkt240_matrix.bin";
const string vector_file = "/scratch/public/sparse_matrix_data/nlpkkt240_vector.bin";

int main(int argc, char** argv)
{
    // Declare forwarders
    extern int mpi(int, char**);
    extern int serial();

    // Usage: program [0|1]
    // (no args) = serial
    // argv[1]: 0 = serial (default), 1 = mpi
    int mode = 0;
    if (argc >= 2) // Check for at least 2 args (program_name + mode)
        mode = atoi(argv[1]); // Read from argv[1]
    
    if (mode == 1) {
        return mpi(argc, argv);
    }

    return serial();
}

// Serial driver implementation (Modified for CSV)
int serial()
{
    CSRMatrix A;
    vector<double> x, y;

    if (!load_matrix(matrix_file, A)) {
        cerr << "Failed to load matrix." << endl;
        return -1;
    }

    if (!load_vector(vector_file, x)) {
        cerr << "Failed to load vector." << endl;
        return -1;
    }

    y.resize(A.n, 0.0);
    auto start = chrono::high_resolution_clock::now();
    
    for (int i = 0; i < A.n; i++) {
        y[i] = 0.0;
        for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; j++)
            y[i] += A.values[j] * x[A.col_idx[j]];
    }
    
    auto end = chrono::high_resolution_clock::now();

    // Use floating-point milliseconds for precision
    chrono::duration<double, milli> duration_ms = end - start;

    // Print CSV format: number_of_processes,computation_time_ms
    cout << "1," << duration_ms.count() << endl;

    return 0;
}

// Build a contiguous row partition that balances by nonzeros. (Unchanged)
static vector<int> build_row_partition_by_nnz(const CSRMatrix& A, int P)
{
    vector<long long> row_nnz(A.n);
    for (int i = 0; i < A.n; ++i) row_nnz[i] = (long long)(A.row_ptr[i+1] - A.row_ptr[i]);

    long long total = 0;
    for (auto v : row_nnz) total += v;
    long long target = (P > 0) ? (total + P - 1) / P : total; // ceil division

    vector<int> starts;
    starts.reserve(P + 1);
    starts.push_back(0);
    long long acc = 0;
    int cur_row = 0;
    for (int r = 1; r < P; ++r) {
        long long goal = (long long)r * target;
        while (cur_row < A.n && acc + row_nnz[cur_row] <= goal) {
            acc += row_nnz[cur_row];
            ++cur_row;
        }
        starts.push_back(cur_row);
    }
    starts.push_back(A.n);

    for (int i = 1; i < (int)starts.size(); ++i) {
        if (starts[i] < starts[i-1]) starts[i] = starts[i-1];
        if (starts[i] > A.n) starts[i] = A.n;
    }
    return starts;
}

// =========================================================================
// MPI IMPLEMENTATION (Modified for CSV)
// =========================================================================
int mpi(int argc, char** argv)
{
    MPI_Init(&argc, &argv); 

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 

    CSRMatrix A;       
    vector<double> x;  
    int n = 0, nnz = 0;
    vector<int> full_row_ptr; 

    // --- Phase 1: Root reads/broadcasts metadata and vector ---
    if (rank == 0) {
        ifstream file(matrix_file, ios::binary);
        if (!file.is_open()) {
            cerr << "Root failed to open matrix file." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file.read(reinterpret_cast<char*>(&n), sizeof(int));
        file.read(reinterpret_cast<char*>(&nnz), sizeof(int));
        
        full_row_ptr.resize(n + 1);
        file.read(reinterpret_cast<char*>(full_row_ptr.data()), (n + 1) * sizeof(int));
        file.close();

        if (!load_vector(vector_file, x)) {
            cerr << "Root failed to load vector." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) full_row_ptr.resize(n + 1);
    if (rank != 0) x.resize(n);
    MPI_Bcast(full_row_ptr.data(), n + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(x.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (n <= 0) {
        if (rank == 0) cerr << "Invalid matrix size." << endl;
        MPI_Finalize();
        return 1;
    }

    // --- Phase 2: Parallel Partitioning ---
    CSRMatrix A_meta;
    A_meta.n = n;
    A_meta.nnz = nnz;
    A_meta.row_ptr = full_row_ptr; 
    vector<int> row_starts = build_row_partition_by_nnz(A_meta, size);
    
    int local_row_start = row_starts[rank];
    int local_row_end = row_starts[rank+1];
    int local_n = local_row_end - local_row_start;
    int global_nnz_offset = full_row_ptr[local_row_start];
    int local_nnz = full_row_ptr[local_row_end] - global_nnz_offset;

    // --- Phase 3: Parallel I/O ---
    long long col_idx_file_start = (long long)sizeof(int) * 2 + (long long)sizeof(int) * (n + 1);
    long long values_file_start = col_idx_file_start + (long long)sizeof(int) * nnz;

    ifstream file(matrix_file, ios::binary);
    if (!file.is_open()) {
        cerr << "Rank " << rank << " failed to open matrix file." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    A.col_idx.resize(local_nnz);
    long long my_col_seek = col_idx_file_start + (long long)sizeof(int) * global_nnz_offset;
    file.seekg(my_col_seek);
    file.read(reinterpret_cast<char*>(A.col_idx.data()), local_nnz * sizeof(int));
    
    A.values.resize(local_nnz);
    long long my_val_seek = values_file_start + (long long)sizeof(double) * global_nnz_offset;
    file.seekg(my_val_seek);
    file.read(reinterpret_cast<char*>(A.values.data()), local_nnz * sizeof(double));
    file.close();

    A.n = local_n;
    A.nnz = local_nnz;
    A.row_ptr.resize(local_n + 1);
    for (int i = 0; i <= local_n; ++i) {
        A.row_ptr[i] = full_row_ptr[local_row_start + i] - global_nnz_offset;
    }

    // --- Phase 4: Compute ---
    vector<double> y_local(local_n, 0.0);
    
    // Barrier to synchronize before computation timer
    MPI_Barrier(MPI_COMM_WORLD);
    double t0_comp = MPI_Wtime(); 

    if (local_n > 0) {
        for (int i = 0; i < local_n; ++i) {
            double sum = 0.0;
            int row_begin = A.row_ptr[i];
            int row_end = A.row_ptr[i+1];
            for (int jj = row_begin; jj < row_end; ++jj) {
                sum += A.values[jj] * x[A.col_idx[jj]];
            }
            y_local[i] = sum;
        }
    }
    
    // Barrier to synchronize after computation
    MPI_Barrier(MPI_COMM_WORLD);
    double t1_comp = MPI_Wtime(); 

    // --- Phase 5: Gather Results ---
    vector<int> recv_counts, recv_displs;
    int my_rows = local_n;
    if (rank == 0) recv_counts.resize(size);
    MPI_Gather(&my_rows, 1, MPI_INT,
               (rank == 0 ? recv_counts.data() : nullptr), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    vector<double> y; 
    if (rank == 0) {
        recv_displs.resize(size, 0);
        for (int r = 1; r < size; ++r) recv_displs[r] = recv_displs[r-1] + recv_counts[r-1];
        y.resize(n, 0.0);
    }
    MPI_Gatherv(y_local.data(), my_rows, MPI_DOUBLE,
                (rank == 0 ? y.data() : nullptr),
                (rank == 0 ? recv_counts.data() : nullptr),
                (rank == 0 ? recv_displs.data() : nullptr), MPI_DOUBLE,
                0, MPI_COMM_WORLD);


    // --- Print Timings ---
    if (rank == 0)
    {
        double comp_ms = (t1_comp - t0_comp) * 1000.0;
        // Print CSV format: number_of_processes,computation_time_ms
        cout << size << "," << comp_ms << endl;
    }

    MPI_Finalize();
    return 0;
}