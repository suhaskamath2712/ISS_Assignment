/*
================================================================================
Sparse Matrix-Vector Multiplication (SpMV) with MPI Parallelization
================================================================================

This program demonstrates both serial and distributed-memory (MPI) parallel sparse
matrix-vector multiplication using the Compressed Sparse Row (CSR) format. It is
designed for large matrices that cannot fit on a single node, and for benchmarking
parallel performance and load balancing.

Key Features:
-------------
* Serial and MPI-based SpMV in a single program, selectable by command-line argument.
* Matrix and vector are loaded from binary files (paths hardcoded for cluster use).
* MPI implementation partitions the matrix rows among processes to balance the number
    of nonzeros per process (for better load balancing).
* Each process receives only its partition of the matrix and computes its local result.
* Results are gathered on the root process for output and timing.

How It Works:
-------------
1. The program is launched with either serial or MPI mode (argv[2]: 0=serial, 1=mpi).
2. In serial mode, the entire matrix and vector are loaded and SpMV is performed.
3. In MPI mode:
     a. The root process loads the matrix/vector and partitions the rows by nonzeros.
     b. The root sends each process its partition (metadata and data) using MPI_Send.
     c. All processes receive the input vector via MPI_Bcast.
     d. Each process computes its local SpMV result.
     e. Results are gathered on the root process using MPI_Gatherv.
     f. The root prints sample results and timing.

Design Notes:
-------------
* The program is robust to any number of MPI processes (including more processes than rows).
* Zero-row processes are handled safely.
* The code is heavily commented for clarity and educational use.
* The matrix/vector file paths are hardcoded for a specific cluster; change as needed.

Usage Examples:
---------------
Serial:   ./main 0 0
MPI:      mpiexec -n 4 ./main 0 1

Arguments:
    argv[1]: (unused, placeholder for future thread count or other options)
    argv[2]: 0 = serial, 1 = MPI

Dependencies:
-------------
* Requires MPI (e.g., OpenMPI, MPICH) for parallel mode.
* Requires a compatible C++17 compiler.
* Matrix/vector binary files must be accessible at the specified paths.

Author: Suhas
Date: 11/11/2025
================================================================================
*/

#include "sparse_matrix.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <mpi.h>

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

    cout << "Serial matrix vector multiplication completed in " << duration.count() << " ms." << endl;
}

// Build a contiguous row partition that balances by nonzeros.
// Returns row_starts of size (P+1) where rows in [row_starts[r], row_starts[r+1]) go to rank r
static vector<int> build_row_partition_by_nnz(const CSRMatrix& A, int P) {
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

int main(int argc, char** argv)
{
    // Declare forwarders
    extern int mpi(int, char**);
    extern int serial();

    // Usage: program <int> [0|1]
    // argv[1]: integer argument (unused)
    // argv[2]: 0 = serial (default), 1 = mpi
    int mode = 0;
    if (argc >= 3) {
        mode = atoi(argv[2]);
    }
    if (mode == 1) return mpi(argc, argv);
    return serial();
}

// Serial driver implementation
int serial()
{
    const string matrix_file = "/scratch/public/sparse_matrix_data/nlpkkt240_matrix.bin";
    const string vector_file = "/scratch/public/sparse_matrix_data/nlpkkt240_vector.bin";

    cout << "=== NLPKKT240 Sparse Matrix-Vector Multiplication Example (Serial) ===" << endl;

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

    print_matrix_info(A);
    y.resize(A.n, 0.0);

    cout << "\nPerforming serial sparse matrix-vector multiplication..." << endl;
    spmv_serial(A, x, y);

    cout << "\nSample results:" << endl;
    cout << "y[0] = " << y[0] << endl;
    cout << "y[" << A.n/2 << "] = " << y[A.n/2] << endl;
    cout << "y[" << A.n-1 << "] = " << y[A.n-1] << endl;

    return 0;
}

// MPI driver implementation
// MPI-based distributed sparse matrix-vector multiplication
// Each MPI process receives a partition of the matrix and computes its local result
int mpi(int argc, char** argv)
{
    // File paths for matrix and vector (must be accessible to all ranks)
    const string matrix_file = "/scratch/public/sparse_matrix_data/nlpkkt240_matrix.bin";
    const string vector_file = "/scratch/public/sparse_matrix_data/nlpkkt240_vector.bin";

    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    int rank = 0, size = 1;
    // Get the rank (process ID) and total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    CSRMatrix A; // Matrix (only fully loaded on root)
    vector<double> x; // Input vector (broadcast to all ranks)

    int n = 0; // Matrix dimension
    vector<int> row_starts; // Row partition boundaries (only on root)

    double t0_total = MPI_Wtime(); // Start total timer

    // Root process loads the matrix and vector, partitions rows
    if (rank == 0) {
        cout << "=== MPI Sparse Matrix-Vector Multiplication (CSR, row-balanced by nnz) ===" << endl;
        if (!load_matrix(matrix_file, A)) {
            cerr << "Failed to load matrix." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (!load_vector(vector_file, x)) {
            cerr << "Failed to load vector." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if ((int)x.size() != A.n) {
            cerr << "Vector length does not match matrix dimension" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        n = A.n;
        // Partition rows among processes for load balancing (by nonzeros)
        row_starts = build_row_partition_by_nnz(A, size);
    }

    // Broadcast matrix dimension to all ranks
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (n <= 0) {
        if (rank == 0) cerr << "Invalid matrix size" << endl;
        MPI_Finalize();
        return 1;
    }

    // Broadcast input vector to all ranks
    if (rank != 0) x.resize(n);
    MPI_Bcast(x.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Each process receives its local matrix partition from root
    int local_row_start = 0, local_row_end = 0;
    int local_n = 0, local_nnz = 0;
    int global_row_offset = 0;
    vector<int> local_row_ptr;
    vector<int> local_col_idx;
    vector<double> local_values;

    if (rank == 0) {
        // Root sends each process its partition (except itself)
        for (int r = 0; r < size; ++r) {
            int rs = row_starts[r];
            int re = row_starts[r+1];
            int rows = re - rs;
            int base = A.row_ptr[rs];
            int lnnz = A.row_ptr[re] - base;

            vector<int> lrp(rows + 1);
            for (int i = 0; i <= rows; ++i) lrp[i] = A.row_ptr[rs + i] - base;

            if (r == 0) {
                // Root keeps its own partition
                local_row_start = rs;
                local_row_end = re;
                local_n = rows;
                local_nnz = lnnz;
                global_row_offset = rs;
                local_row_ptr = std::move(lrp);
                local_col_idx.assign(A.col_idx.begin() + base, A.col_idx.begin() + base + lnnz);
                local_values.assign(A.values.begin() + base, A.values.begin() + base + lnnz);
            } else {
                // Send partition metadata and data to rank r
                MPI_Send(&rs, 1, MPI_INT, r, 100, MPI_COMM_WORLD); // Row start
                MPI_Send(&re, 1, MPI_INT, r, 101, MPI_COMM_WORLD); // Row end
                MPI_Send(&rows, 1, MPI_INT, r, 102, MPI_COMM_WORLD); // Number of rows
                MPI_Send(&lnnz, 1, MPI_INT, r, 103, MPI_COMM_WORLD); // Number of nonzeros
                MPI_Send(lrp.data(), rows + 1, MPI_INT, r, 104, MPI_COMM_WORLD); // Row pointer
                MPI_Send(A.col_idx.data() + base, lnnz, MPI_INT, r, 105, MPI_COMM_WORLD); // Column indices
                MPI_Send(A.values.data() + base, lnnz, MPI_DOUBLE, r, 106, MPI_COMM_WORLD); // Values
            }
        }
    } else {
        // Non-root ranks receive their partition from root
        MPI_Recv(&local_row_start, 1, MPI_INT, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_row_end, 1, MPI_INT, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_n, 1, MPI_INT, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_nnz, 1, MPI_INT, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        global_row_offset = local_row_start;

        local_row_ptr.resize(local_n + 1);
        local_col_idx.resize(local_nnz);
        local_values.resize(local_nnz);
        MPI_Recv(local_row_ptr.data(), local_n + 1, MPI_INT, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Row pointer
        MPI_Recv(local_col_idx.data(), local_nnz, MPI_INT, 0, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Column indices
        MPI_Recv(local_values.data(), local_nnz, MPI_DOUBLE, 0, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Values
    }

    vector<double> y_local(local_n, 0.0);
    double t0_comp = MPI_Wtime();
    if (local_n > 0) {
        for (int i = 0; i < local_n; ++i) {
            double sum = 0.0;
            int row_begin = local_row_ptr[i];
            int row_end = local_row_ptr[i+1];
            for (int jj = row_begin; jj < row_end; ++jj) {
                sum += local_values[jj] * x[local_col_idx[jj]];
            }
            y_local[i] = sum;
        }
    }
    double t1_comp = MPI_Wtime();

    vector<int> recv_counts, recv_displs;
    int my_rows = local_n;
    if (rank == 0) recv_counts.resize(size);
    MPI_Gather(&my_rows, 1, MPI_INT,
               (rank == 0 ? recv_counts.data() : nullptr), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    vector<double> y; // full result on root
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

    double t1_total = MPI_Wtime();

    if (rank == 0)
    {
        cout << "\nSample results:" << endl;
        if (!y.empty()) {
            cout << "y[0] = " << y[0] << endl;
            cout << "y[" << n/2 << "] = " << y[n/2] << endl;
            cout << "y[" << n-1 << "] = " << y[n-1] << endl;
        }
        double total_ms = (t1_total - t0_total) * 1000.0;
        cout << "Total elapsed (wall): " << total_ms << " ms across " << size << " rank(s)." << endl;
        long long total_rows = 0; for (int c : recv_counts) total_rows += c;
        cout << "Distributed rows: " << total_rows << " (matrix rows: " << n << ")" << endl;
        if ((int)recv_counts.size() != size) cout << "[Debug] recv_counts size mismatch" << endl;
    }

    MPI_Finalize();
    return 0;
}
