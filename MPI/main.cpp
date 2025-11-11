/*
================================================================================
Sparse Matrix-Vector Multiplication (SpMV) with MPI Parallelization
================================================================================
Arguments:
    argv[1]: 0 = serial (default), 1 = mpi
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

// Serial driver implementation
int serial()
{
    CSRMatrix A;
    vector<double> x, y;

    if (!load_matrix(matrix_file, A)) {
        cerr << "Failed to load matrix." << endl;
        cerr << "Make sure you have read access to /scratch/public/sparse_matrix_data/" << endl;
        return -1;
    }

    if (!load_vector(vector_file, x)) {
        cerr << "Failed to load vector." << endl;
        return -1;
    }

    y.resize(A.n, 0.0);

    cout << "Mode: Serial" << endl;
    auto start = chrono::high_resolution_clock::now();
    
    for (int i = 0; i < A.n; i++) {
        y[i] = 0.0;
        for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; j++)
            y[i] += A.values[j] * x[A.col_idx[j]];
    }
    
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << "Serial matrix vector multiplication completed in " << duration.count() << " ms." << endl;

    return 0;
}

// Build a contiguous row partition that balances by nonzeros.
// Returns row_starts of size (P+1) where rows in [row_starts[r], row_starts[r+1]) go to rank r
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

// MPI-based distributed sparse matrix-vector multiplication
int mpi(int argc, char** argv)
{
    MPI_Init(&argc, &argv); 

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 

    if (rank == 0) {
        cout << "Mode: MPI (running on " << size << " processes)" << endl;
    }

    CSRMatrix A; 
    vector<double> x; 
    int n = 0; 
    vector<int> row_starts; 

    double t0_total = MPI_Wtime(); 

    struct RankData {
        int rs, re, rows, lnnz;
        vector<int> lrp;
    };
    
    vector<RankData> all_rank_data;
    vector<MPI_Request> send_requests;

    if (rank == 0)
    {
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
        row_starts = build_row_partition_by_nnz(A, size); 
        all_rank_data.resize(size);

        for (int r = 0; r < size; ++r) {
            auto& data = all_rank_data[r]; 
            data.rs = row_starts[r];
            data.re = row_starts[r+1];
            data.rows = data.re - data.rs;
            int base = A.row_ptr[data.rs];
            data.lnnz = A.row_ptr[data.re] - base;

            data.lrp.resize(data.rows + 1);
            for (int i = 0; i <= data.rows; ++i) {
                data.lrp[i] = A.row_ptr[data.rs + i] - base;
            }
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (n <= 0)
    {
        if (rank == 0) cerr << "Invalid matrix size" << endl;
        MPI_Finalize();
        return 1;
    }

    if (rank != 0) x.resize(n);
    MPI_Bcast(x.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int local_row_start = 0, local_row_end = 0;
    int local_n = 0, local_nnz = 0;
    vector<int> local_row_ptr;
    vector<int> local_col_idx;
    vector<double> local_values;

    if (rank == 0) {
        for (int r = 1; r < size; ++r) { 
            auto& data = all_rank_data[r]; 
            int base = A.row_ptr[data.rs];
            
            size_t start_req_idx = send_requests.size();
            send_requests.resize(start_req_idx + 7);

            MPI_Isend(&data.rs, 1, MPI_INT, r, 100, MPI_COMM_WORLD, &send_requests[start_req_idx + 0]);
            MPI_Isend(&data.re, 1, MPI_INT, r, 101, MPI_COMM_WORLD, &send_requests[start_req_idx + 1]);
            MPI_Isend(&data.rows, 1, MPI_INT, r, 102, MPI_COMM_WORLD, &send_requests[start_req_idx + 2]);
            MPI_Isend(&data.lnnz, 1, MPI_INT, r, 103, MPI_COMM_WORLD, &send_requests[start_req_idx + 3]);
            MPI_Isend(data.lrp.data(), data.rows + 1, MPI_INT, r, 104, MPI_COMM_WORLD, &send_requests[start_req_idx + 4]);
            MPI_Isend(A.col_idx.data() + base, data.lnnz, MPI_INT, r, 105, MPI_COMM_WORLD, &send_requests[start_req_idx + 5]);
            MPI_Isend(A.values.data() + base, data.lnnz, MPI_DOUBLE, r, 106, MPI_COMM_WORLD, &send_requests[start_req_idx + 6]);
        }

        auto& data = all_rank_data[0];
        int base = A.row_ptr[data.rs];
        local_row_start = data.rs;
        local_row_end = data.re;
        local_n = data.rows;
        local_nnz = data.lnnz;
        // global_row_offset = data.rs;  <-- REMOVED
        local_row_ptr = std::move(data.lrp); 
        local_col_idx.assign(A.col_idx.begin() + base, A.col_idx.begin() + base + data.lnnz);
        local_values.assign(A.values.begin() + base, A.values.begin() + base + data.lnnz);
    }
    else
    {
        MPI_Recv(&local_row_start, 1, MPI_INT, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_row_end, 1, MPI_INT, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_n, 1, MPI_INT, 0, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_nnz, 1, MPI_INT, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // global_row_offset = local_row_start; <-- REMOVED

        local_row_ptr.resize(local_n + 1);
        local_col_idx.resize(local_nnz);
        local_values.resize(local_nnz);
        MPI_Recv(local_row_ptr.data(), local_n + 1, MPI_INT, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_col_idx.data(), local_nnz, MPI_INT, 0, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_values.data(), local_nnz, MPI_DOUBLE, 0, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if (rank == 0 && !send_requests.empty())
        MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

    vector<double> y_local(local_n, 0.0);
    double t0_comp = MPI_Wtime(); // <-- Timer start
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
    double t1_comp = MPI_Wtime(); // <-- Timer end

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

    double t1_total = MPI_Wtime();

    if (rank == 0)
    {
        double total_ms = (t1_total - t0_total) * 1000.0;
        cout << "Total elapsed: " << total_ms << " ms." << endl;
        // --- FIX: Print the computation time ---
        cout << "  (Computation time: " << (t1_comp - t0_comp) * 1000.0 << " ms)" << endl;
    }

    MPI_Finalize();
    return 0;
}