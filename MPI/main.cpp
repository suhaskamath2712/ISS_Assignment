/*
================================================================================
Sparse Matrix-Vector Multiplication (SpMV) with MPI Parallelization
================================================================================
This version outputs per-process computation times and total run times.
CSV Format:
#,number_of_processes,rank_of_process,computation_time_ms
$,number_of_processes,total_execution_time_ms
================================================================================
*/

#include "sparse_matrix.h"
#include <iostream>
#include <vector>
#include <chrono> // For serial timing
#include <algorithm>
#include <numeric>
#include <fstream> // Needed for file I/O in parallel
#include <mpi.h>

//Preprocessor macro for ceiling division
#define CEIL_DIV(a, b) (((a) + (b) - 1) / (b))

//Preprocessor macro to add commas for CSV
#define COMMA ","

using namespace std;

//File paths for matrix and vector 
const string matrix_file = "/scratch/public/sparse_matrix_data/nlpkkt240_matrix.bin";
const string vector_file = "/scratch/public/sparse_matrix_data/nlpkkt240_vector.bin";

// --- Forward Declarations ---
int serial();
int mpi(int argc, char** argv);
int read_metadata(int& n, int& nnz, vector<int>& full_row_ptr, vector<double>& x);
int read_local_data(int rank, int n, int nnz, int local_n, int global_nnz_offset, int local_nnz, CSRMatrix& A);
static vector<int> build_row_partition_by_nnz(int n, int nnz, const vector<int>& full_row_ptr, int P);


// --- Main Function ---
int main(int argc, char** argv)
{
    // Usage: program [0|1]
    // (no args) = serial
    // argv[1]: 0 = serial (default), 1 = mpi
    int mode = 0;
    if (argc >= 2)
        mode = atoi(argv[1]);
    
    if (mode == 1) {
        return mpi(argc, argv);
    }

    return serial();
}

// --- Serial Implementation ---
int serial()
{
    // Timer for total execution
    auto start_total = chrono::high_resolution_clock::now();

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
    
    // Timer for computation
    auto start_comp = chrono::high_resolution_clock::now();
    
    for (int i = 0; i < A.n; i++) {
        y[i] = 0.0;
        for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; j++)
            y[i] += A.values[j] * x[A.col_idx[j]];
    }
    
    auto end_comp = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> comp_ms = end_comp - start_comp;

    // Print computation time
    cout << "#" << COMMA << 1 << COMMA << 0 << COMMA << comp_ms.count() << endl;

    // Stop total timer and print total time
    auto end_total = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> total_ms = end_total - start_total;

    // Print total time
    cout << "$" << COMMA << 1 << COMMA << total_ms.count() << endl;

    return 0;
}

// --- MPI Implementation ---
int mpi(int argc, char** argv)
{
    MPI_Init(&argc, &argv); 

    // Start total timer
    double t0_total = MPI_Wtime();

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 

    CSRMatrix A;       
    vector<double> x;  
    int n = 0, nnz = 0;
    vector<int> full_row_ptr; 

    // --- Phase 1: Root reads/broadcasts metadata and vector ---
    if (rank == 0) {
        if (read_metadata(n, nnz, full_row_ptr, x) != 0) {
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
    vector<int> row_starts = build_row_partition_by_nnz(n, nnz, full_row_ptr, size);
    
    int local_row_start = row_starts[rank];
    int local_row_end = row_starts[rank+1];
    int local_n = local_row_end - local_row_start;
    int global_nnz_offset = full_row_ptr[local_row_start];
    int local_nnz = full_row_ptr[local_row_end] - global_nnz_offset;

    // --- Phase 3: Parallel I/O ---
    if (read_local_data(rank, n, nnz, local_n, global_nnz_offset, local_nnz, A) != 0) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Fix up the local row_ptr (0-based)
    A.row_ptr.resize(local_n + 1);
    for (int i = 0; i <= local_n; ++i) {
        A.row_ptr[i] = full_row_ptr[local_row_start + i] - global_nnz_offset;
    }

    // --- Phase 4: Compute ---
    vector<double> y_local(local_n, 0.0);
    
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
    
    // All processes print their computation time
    double comp_ms = (t1_comp - t0_comp) * 1000.0;
    cout << "#" << COMMA << size << COMMA << rank << COMMA << comp_ms << endl;

    // Stop total timer
    double t1_total = MPI_Wtime();
    
    // Only Rank 0 prints the total execution time

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        double total_ms = (t1_total - t0_total) * 1000.0;
        // Print total time
        cout << "$" << COMMA << size << COMMA << total_ms << endl;
    }

    MPI_Finalize();
    return 0;
}


// --- Helper Functions ---

// Build a contiguous row partition that balances by nonzeros.
static vector<int> build_row_partition_by_nnz(int n, int nnz, const vector<int>& full_row_ptr, int P)
{
    vector<long long> row_nnz(n);
    for (int i = 0; i < n; ++i) row_nnz[i] = (long long)(full_row_ptr[i+1] - full_row_ptr[i]);

    long long total = 0;
    for (auto v : row_nnz) total += v;
    long long target = (P > 0) ? (total + P - 1) / P : 0; // Use CEIL_DIV logic

    vector<int> starts;
    starts.reserve(P + 1);
    starts.push_back(0);
    long long acc = 0;
    int cur_row = 0;
    for (int r = 1; r < P; ++r) {
        long long goal = (long long)r * target;
        while (cur_row < n && acc + row_nnz[cur_row] <= goal) {
            acc += row_nnz[cur_row];
            ++cur_row;
        }
        starts.push_back(cur_row);
    }
    starts.push_back(n);

    for (int i = 1; i < (int)starts.size(); ++i) {
        if (starts[i] < starts[i-1]) starts[i] = starts[i-1];
        if (starts[i] > n) starts[i] = n;
    }
    return starts;
}

// Function for Rank 0 to read metadata
int read_metadata(int& n, int& nnz, vector<int>& full_row_ptr, vector<double>& x)
{
    ifstream file(matrix_file, ios::binary);
    if (!file.is_open()) {
        cerr << "Root failed to open matrix file." << endl;
        return -1;
    }
    file.read(reinterpret_cast<char*>(&n), sizeof(int));
    file.read(reinterpret_cast<char*>(&nnz), sizeof(int));
    
    full_row_ptr.resize(n + 1);
    file.read(reinterpret_cast<char*>(full_row_ptr.data()), (n + 1) * sizeof(int));
    file.close();

    if (!load_vector(vector_file, x)) {
        cerr << "Root failed to load vector." << endl;
        return -1;
    }
    return 0;
}

// Function for all ranks to read their local data in parallel
int read_local_data(int rank, int n, int nnz, int local_n, int global_nnz_offset, int local_nnz, CSRMatrix& A)
{
    // Calculate file offsets
    long long col_idx_file_start = (long long)sizeof(int) * 2 + (long long)sizeof(int) * (n + 1);
    long long values_file_start = col_idx_file_start + (long long)sizeof(int) * nnz;

    ifstream file(matrix_file, ios::binary);
    if (!file.is_open()) {
        cerr << "Rank " << rank << " failed to open matrix file." << endl;
        return -1;
    }

    // Read local col_idx
    A.col_idx.resize(local_nnz);
    long long my_col_seek = col_idx_file_start + (long long)sizeof(int) * global_nnz_offset;
    file.seekg(my_col_seek);
    file.read(reinterpret_cast<char*>(A.col_idx.data()), local_nnz * sizeof(int));
    
    // Read local values
    A.values.resize(local_nnz);
    long long my_val_seek = values_file_start + (long long)sizeof(double) * global_nnz_offset;
    file.seekg(my_val_seek);
    file.read(reinterpret_cast<char*>(A.values.data()), local_nnz * sizeof(double));
    file.close();

    A.n = local_n;
    A.nnz = local_nnz;
    return 0;
}